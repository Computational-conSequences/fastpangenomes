#!/usr/bin/env perl

#########################################################################
#                                                                       #
#       Author: Gabo Moreno-Hagelsieb                                   #
#                                                                       #
#########################################################################
use Getopt::Long;
use strict;
# to make temporary files/directories
use File::Temp qw( tempfile tempdir );
use sigtrap qw(handler signalHandler normal-signals);

### this is to be able to add an option for the
### clustering method.
# diana is a divisive
my @aggMethods = qw(
                         average
                         single
                         complete
                         mcquitty
                         ward
                 );

my $clustMethodsMatch = join("|","diana",@aggMethods);
my $aggMethods        = join("|",@aggMethods);

my @metrics = qw(
                    sign
                    gss
                    ani
                    mash
                    mdashing
                    dashing
                    jaccard
            );
my $metricMatch = join("|",@metrics);

my %metrics = (
    "sign"     => qq(DNA signatures),
    "gss"      => qq(Genomic Similarity Scores),
    "ani"      => qq(Average Nucleotide Identity),
    "mash"     => qq(mash's MinHash distance),
    "mdashing" => qq(dashing's MinHash distance),
    "dashing"  => qq(dashing's Jaccard distance),
    "jaccard"  => qq(domain matches Jaccard distance),
);

######### defaults
my $calcFile     = '';
my $listFile     = 'all';
my $outputFolder = 'Clusters';
my $clustMethod  = 'diana';
my $defBasis     = 'sign';
my $inBasis      = '';
my $cuts         = 0;

my $ownName = $0;
$ownName =~ s{.*/}{};
my $helpMsg
    = qq(usage:\n)
    . qq(   $ownName [options]\n)
    . qq(\n)
    . qq(options:\n)
    . qq(   -i table with genome distances/similarities, such as\n)
    . qq(       ANI, mash, dashing, mash-dashing, DNA signatures,\n)
    . qq(       GSS, domain abundances, required\n)
    . qq(   -b metric type [$metricMatch],\n)
    . qq(       default sign, or determined from the file's name:\n)
    . qq(       sign: DNA signature table (each row is a DNA signature)\n)
    . qq(       gss:  Genomic Similarity Score (fraction), pairwise table\n)
    . qq(       ani:  Average Nucleotide Identity (percent), pairwise table\n)
    . qq(       mash: mash distances, pairwise table\n)
    . qq(       mdashing: dashing's mash distances matrix\n)
    . qq(       dashing:  dashing's Jaccard index matrix\n)
    . qq(       jaccard:  domain abundance matrix\n)
    . qq(   -l file with list of genomes to cluster,\n)
    . qq(       default all genomes [all]\n)
    . qq(   -o output folder, default: Clusters\n)
    . qq(   -m clustering method: divisive [diana], or\n)
    . qq(       agglomerative [$aggMethods],\n)
    . qq(       default: $clustMethod\n)
    . qq(   -c cut into genome groups, 0 (none), 1 (all),\n)
    . qq(       specific (e.g. 0.03), default: $cuts\n)
    . qq(\n)
    . qq(requirements:\n)
    . qq(    R from r-project.org\n)
    . qq(    R packages: cluster, fastcluster, MCMCpack, ape and reshape2\n\n)
    ;


my $options = GetOptions(
    "i=s" => \$calcFile,
    "l=s" => \$listFile,
    "o=s" => \$outputFolder,
    "m=s" => \$clustMethod,
    "b=s" => \$inBasis,
    "c=s" => \$cuts,
) or die $helpMsg;

if ( !$calcFile ) {
    die $helpMsg;
}

### file is either signatures, GSS, or ani
my $basis
    = $inBasis  =~ m{$metricMatch}i                        ? lc($&)
    : $calcFile =~ m{dashing}i && $calcFile =~ m{mash}i    ? "mdashing"
    : $calcFile =~ m{dashing}i && $calcFile =~ m{jaccard}i ? "dashing"
    : $calcFile =~ m{mash}i                                ? "mash"
    : $calcFile =~ m{ani}i                                 ? "ani"
    : $calcFile =~ m{gss}i                                 ? "gss"
    : $calcFile =~ m{pfam|cdd|vfdb}i                       ? "jaccard"
    : 'sign';

### test that the method is spelled correctly:

my $procedure = '';
if( $clustMethod =~ m{^($clustMethodsMatch)$}i ) {
    $clustMethod = lc($1);
    $clustMethod =~ s{\.d}{.D};
    my $addprocedure
        = $clustMethod eq "ward"
        ? qq(fastcluster::hclust(dx,method="ward.D2"))
        : $clustMethod eq "diana"
        ? qq(diana(dx))
        : qq(fastcluster::hclust(dx,method="$clustMethod"));
    $procedure .= $addprocedure;
    print "clustering genomes by ",$metrics{"$basis"}," ($basis),\n";
    print "   using $clustMethod clustering method\n";
}
else {
    die "   the clustering method [$clustMethod] does not exist\n"
        . "   try any of [$clustMethodsMatch] (default: diana)\n\n";
}

######################################################################
############# Start working ##########################################
######################################################################
#        library("MCMCpack")
my $Rhead = << 'RHEAD100';
### R script to cluster and cut genomes
suppressWarnings(
    suppressPackageStartupMessages({
        library("cluster")
        library("MCMCpack")
        library("ape")
        library("reshape2")
    })
)
### this is a function to cut clusters at a given threshold.
### Inputs are dendrogram (tree) and threshold of interest.
prune<-function(dendrogram,threshold) {
    g<-cutree(dendrogram,h=threshold)
    groupnumber<-length(table(g))
    engorda<-paste('%0',nchar(groupnumber),"d",sep="")
    b<-noquote(paste("Group",sprintf(engorda,(1:groupnumber)),sep = "-"))
    namelists=split(names(g),g)
    data.frame(noquote(I(unlist(lapply(namelists,paste,collapse=",")))),
               row.names=noquote(b))
}
### done with the function
RHEAD100

my $cwd = qx(pwd);
chomp($cwd);
my $calcFile
    = $calcFile =~ m{(^/)} ? $calcFile
    : $cwd . "/" . $calcFile;
unless( -f "$calcFile" ) {
    die "I cannot find\n$calcFile\n\n";
}

### test the list file if other than all
my %good = ();
unless( $listFile eq "all" ) {
    print "      learning IDs\n";
    if( -f "$listFile" ) {
        open( my $LIST,"<","$listFile" );
        while(<$LIST>) {
            chomp;
            s{^\s+}{};
            s{\s+$}{};
            my $genome = $_;
            if( length($genome) > 2 ) {
                $good{"$genome"}++;
            }
        }
        my $good = keys %good;
        if( $good == 0 ) {
            die "$listFile seems to be empty\n\n";
        }
        else {
            print "will cluster $good genomes\n";
        }
    }
    else {
        die "$listFile does not exist\n\n";
    }
}

### check cuts
my $groups
    = $cuts > 0 && $cuts < 1 ? $cuts
    : $cuts > 60 && $cuts < 100 && ( $inBasis eq "ani" ) ? $cuts
    : $cuts == 1 ? "1:10"
    : "none";
if( $groups eq "1:10" ) {
    print "   cutting at a range of thresholds\n";
}
else {
    print "   cutting clusters at a distance of $groups\n";
}

my $tmpDir = tempdir("/tmp/$ownName.XXXXXXXXXXXX");
print "   working in temp folder:\n    $tmpDir\n";

my $listroot = $listFile;
$listroot =~ s{.+\/}{};
$listroot =~ s{\.\w+$}{};
my $calcroot = $calcFile;
$calcroot =~ s{.+\/}{};
$calcroot =~ s{\.\S+$}{};
my $root = $listroot eq "all" ? $calcroot : join("-",$calcroot,$listroot);
print "   root name is $root\n";

### now save distances into file to be read by R

if( $basis eq "mash" ) {
    $calcFile = fixMash("$calcFile");
}
elsif( $basis =~ m{dashing} ) {
    $calcFile = fixDashing("$calcFile");
}
my $readFile
    = $listFile eq "all" ? $calcFile : "$root.subset.bz2";
unless( $readFile eq $calcFile ) {
    print "   extracting subset\n";
    my ($open,$cfile) = figureCompression("$calcFile");
    my $printed = 0;
    open( my $INFILE,"-|","$open $cfile" );
    open( my $SUBFILE,"|-","bzip2 -9 > $tmpDir/$readFile" );
    if( $basis eq "sign"
        || $basis eq "jaccard" ) {
        while(<$INFILE>) {
            my($genome) = split(/\t/,$_);
            if( $. == 1 or $good{"$genome"} > 0 ) {
                print {$SUBFILE} $_;
                $printed++;
            }
        }
    }
    elsif( $basis =~ m{dashing} ) {
        my $keep = '';
        my @keep = ();
      DASHROW:
        while(<$INFILE>) {
            my @kept = ();
            if( $. == 1 ) {
                my @genomes = split;
                for my $i ( 0 .. $#genomes ) {
                    if( $good{"$genomes[$i]"} > 0 ){
                        push(@kept,$genomes[$i]);
                        push(@keep,$i);
                    }
                }
            }
            else {
                my( $genome,@dists ) = split;
                if( $good{"$genome"} > 0 ) {
                    push(@kept,$genome,@dists[@keep]);
                    $printed++;
                }
                else {
                    next DASHROW;
                }
            }
            print {$SUBFILE} join("\t",@kept),"\n";
        }
    }
    else {
        while(<$INFILE>) {
            my($genome1,$genome2) = split(/\t/,$_);
            if( $. == 1 ) {
                print {$SUBFILE} $_;
            }
            elsif(
                exists($good{"$genome1"}) && exists($good{"$genome2"})
            ) {
                print {$SUBFILE} $_;
                $printed++;
            }
        }
    }
    close($INFILE);
    close($SUBFILE);
    if( $printed < 2 ) {
        print "The genomes in $listFile are not in\n$calcFile\n";
        signalHandler();
    }
}

print "   clustering:\n";

my $clustFile = "$root.$clustMethod";
if( $basis eq "sign" ) {
    clustSign("$clustMethod","$readFile","$clustFile","$groups");
    runR("$clustFile");
}
elsif( $basis eq "jaccard" ) {
    for my $weight ( 'weighted','binary' ) {
        my $jaccFile = $root . "-" . $weight . "." . $clustMethod;
        clustJaccard("$clustMethod","$readFile","$jaccFile",
                     "$groups","$weight");
        runR("$jaccFile");
    }
}
elsif( $basis eq "ani" ) {
    clustANI("$clustMethod","$readFile","$clustFile","$groups");
    runR("$clustFile");
}
elsif( $basis eq "mash" ) {
    clustMash("$clustMethod","$readFile","$clustFile","$groups");
    runR("$clustFile");
}
elsif( $basis eq "mash" ) {
    clustMash("$clustMethod","$readFile","$clustFile","$groups");
    runR("$clustFile");
}
elsif( $basis =~ m{dashing} ) {
    clustDashing("$clustMethod","$readFile","$clustFile","$groups");
    runR("$clustFile");
}
else { ### this is GSS
    for my $level ( "a","b","c" ) {
        my $gssFile = "GSS" . $level . "." . $clustMethod;
        clustGSS("$clustMethod","$readFile","$gssFile",
                 "$groups","$level");
        runR("$gssFile");
    }
}

unless( -d "$outputFolder" ) {
    mkdir("$outputFolder");
}
system("rsync -aqz $tmpDir/* $outputFolder 2>/dev/null");
print "   The cluster has been saved in Newick format:\n"
    . "      $outputFolder/$clustFile.tree\n";

print "\tcleaning up ...";
system "rm -r $tmpDir";
print  "\n        $ownName done!\n\n";

############# sub routines #############
sub signalHandler {
    my $message = $_[0];
    if( length("$message") > 0 ) {
        print $message,"\n";
    }
    if( -d "$tmpDir" ) {
        print "\n\tcleaning up ...\n";
        system "rm -r $tmpDir";
        die  "    done!\n\n";
    }
    else {
        print "\n\ttemp files cleared out\n";
        die  "    done!\n\n";
    }
}

#### Cluster by DNA signature
sub clustSign {
    my($clustMethod,$subFile,$outFile,$cutGroups) = @_;
    my $rFile = "$tmpDir/$root.$clustMethod.R";
    print "      producing R file:\n"
        . "$rFile\n";
    open( my $CLUSTSIGN,">","$rFile" );
    print {$CLUSTSIGN} << "CLUSTSIGN";
$Rhead

############## import the distance data
signaturetbl<-read.table("$subFile", sep="\\t")

############## calculate the distance matrix
dx <- dist(signaturetbl,method="manhattan")
### Karlin's delta is normalized against number of nucleotides
nucl <- ncol(signaturetbl)
dx = dx/nucl

#### the output file contains the distance matrix
xmat = as.matrix(dx)
write.table(xmat,
        bzfile("$root.dist.matrix.bz2",compression = 9),
        row.names=T,col.names=T,quote=F)

#### the output file contains the distance matrix as a table
xmat[upper.tri(xmat)] <- NA
x = melt(xmat, na.rm=T)
write.table(x, bzfile("$root.dist.tbl.bz2",compression = 9),
               quote=F,row.names=F,col.names=F,sep="\\t")

## create a dendrogram using the agglomerative (bottom up, agnes)
## or divisive (top-down, diana) and the dx dissimilarity matrix
gnmCluster<-$procedure
gnmDend<-as.hclust(gnmCluster)
saveRDS(gnmDend,"${outFile}-cluster.Rds")

## Output the tree in the Newick format
phy <- as.phylo(as.hclust(gnmCluster))
write.tree(phy, file="$outFile.tree")

## calculate the divisive/agglomerative coefficients
coeff.gnm<-coef(gnmCluster)
coeff.gnm
coeff.print<-formatC(coeff.gnm,digits=3,format="f")
write(paste("$clustMethod", coeff.print),
      "$root.coeff",sep = "\\t", append=T)
CLUSTSIGN
    if( $cutGroups eq "1:10" ) {
        print {$CLUSTSIGN} << "APPENDS";
for ( i in $cutGroups ){
    thrd<-i * 0.01
    thrname = formatC(thrd,digits=2,format="f")
    groups<-prune(gnmDend,thrd)
    write.table(groups,
                file=paste(paste("$outFile",thrname,sep="-"),"groups",sep="."),
                row.names=T,col.names=F,quote=F,append=F)
}
APPENDS
    }
    elsif( $cutGroups ne "none" ) {
        print {$CLUSTSIGN} << "APPENDS1";
thrd = $cutGroups
thrname = formatC(thrd,digits=3,format="f")
groups<-prune(gnmDend,thrd)
write.table(groups,
            file=paste(paste("$outFile",thrname,sep="-"),"groups",sep="."),
            row.names=T,col.names=F,quote=F,append=F)
APPENDS1
    }
    close($CLUSTSIGN);
}

#### Cluster by Genomic Similarity Score (or some fraction thing)
sub clustGSS {
    my($clustMethod,$subFile,$outFile,$cutGroups,$level) = @_;
    my $column
        = $level eq "a" ? 3
        : $level eq "b" ? 4
        : $level eq "c" ? 5
        : '3';
    my $rFile = "$tmpDir/$outFile.R";
    print "      producing R file:\n"
        . "$rFile\n";
    open( my $CLUSTGSS,">","$rFile" );
    print {$CLUSTGSS} << "CLUSTGSS";
$Rhead

############## import GSS data
gsstbl<-read.table("$subFile",header=T,sep="\\t")

##### obtain a vector of genome names to be able to name the
##### items in redundant group later
number.genome<-length(levels(as.factor(gsstbl[,2])))
genome.names<-gsstbl[1:number.genome,2]

##### define the distance matrix of using 1-gsstbl
##### GSS (GSSa = 3, GSSb = 4, GSSc = 5)
dx<-xpnd(gsstbl[,$column])
dx<-1-dx
rownames(dx)<-genome.names
colnames(dx)<-genome.names

##### the output file contains the distance matrix
xmat = as.matrix(dx)
write.table(xmat,
        bzfile("$root.dist.matrix.bz2",compression = 9),
        row.names=T,col.names=T,quote=F)

#### the output file contains the distance matrix as a table
xmat[upper.tri(xmat)] <- NA
x = melt(xmat, na.rm=T)
write.table(x, bzfile("$root.dist.tbl.bz2",compression = 9),
               quote=F,row.names=F,col.names=F,sep="\\t")

##### create a dendrogram using the agglomerative (bottom up, agnes)
##### or divisive (top-down, diana) and the dx dissimilarity matrix
dx<-as.dist(dx)
gnmCluster<-$procedure
gnmDend<-as.hclust(gnmCluster)
saveRDS(gnmDend,"${outFile}-cluster.Rds")

##### Output the tree in the Newick format
phy <- as.phylo(as.hclust(gnmCluster))
write.tree(phy, file="$outFile.tree")

##### calculate the divisive/agglomerative coefficients
coeff.gnm<-coef(gnmCluster)
coeff.gnm
coeff.print<-formatC(coeff.gnm,digits=3,format="f")
write(paste("$clustMethod", coeff.print),
      "$outFile.coeff",sep = "\\t", append=T)
CLUSTGSS
    if( $cutGroups eq "1:10" ) {
        print {$CLUSTGSS} << "APPENDG";
for ( i in $cutGroups ){
    thrd<-i * 0.05
    uthrname = 1 - thrd
    thrname = formatC(uthrname,digits=2,format="f")
    groups<-prune(gnmDend,thrd)
    write.table(groups,
                file=paste(paste("$outFile",thrname,sep="-"),"groups",sep="."),
                row.names=T,col.names=F,quote=F,append=F)
}
APPENDG
    }
    elsif( $cutGroups ne "none" ) {
        print {$CLUSTGSS} << "APPENDG1";
thrd = $cutGroups
uthrname = 1 - thrd
thrname = formatC(uthrname,digits=2,format="f")
groups<-prune(gnmDend,thrd)
write.table(groups,
            file=paste(paste("$outFile",thrname,sep="-"),"groups",sep="."),
            row.names=T,col.names=F,quote=F,append=F)
APPENDG1
    }
    close($CLUSTGSS);
}

##### cluster by ANI (or some percent thing)
sub clustANI {
    my($clustMethod,$subFile,$outFile,$cutGroups) = @_;
    my $rFile = "$tmpDir/$outFile.R";
    print "      producing R file:\n"
        . "$rFile\n";
    open( my $CLUSTANI,">","$rFile" );
    print {$CLUSTANI} << "CLUSTANI";
$Rhead

############## import ANI data
anitbl<-read.table("$subFile",header=T,sep="\\t")

##### obtain a vector of genome names to be able to name the
##### items in redundant group later
number.genome<-length(levels(as.factor(anitbl[,2])))
genome.names<-anitbl[1:number.genome,2]

##### define the distance matrix of using 1-anitbl
dx <- xpnd(anitbl[,ncol(anitbl)])
dx <- (100 - dx)/100
rownames(dx)<-genome.names
colnames(dx)<-genome.names

##### the output file contains the distance matrix
xmat = as.matrix(dx)
write.table(xmat,
        bzfile("$root.dist.matrix.bz2",compression = 9),
        row.names=T,col.names=T,quote=F)

#### the output file contains the distance matrix as a table
xmat[upper.tri(xmat)] <- NA
x = melt(xmat, na.rm=T)
write.table(x, bzfile("$root.dist.tbl.bz2",compression = 9),
               quote=F,row.names=F,col.names=F,sep="\\t")

##### create a dendrogram using the agglomerative (bottom up, agnes)
##### or divisive (top-down, diana) and the dx dissimilarity matrix
dx<-as.dist(dx)
gnmCluster<-$procedure
gnmDend<-as.hclust(gnmCluster)
saveRDS(gnmDend,"${outFile}-cluster.Rds")

##### Output the tree in the Newick format
phy <- as.phylo(as.hclust(gnmCluster))
write.tree(phy, file="$outFile.tree")

##### calculate the divisive/agglomerative coefficients
coeff.gnm<-coef(gnmCluster)
coeff.gnm
coeff.print<-formatC(coeff.gnm,digits=3,format="f")
write(paste("$clustMethod", coeff.print),
      "$outFile.coeff",sep = "\\t", append=T)
CLUSTANI
    if( $cutGroups eq "1:10" ) {
        print {$CLUSTANI} << "APPENDANI";
for ( i in 1:10 ){
    thrd<-i * 0.5
    uthrname = 100 - thrd
    thrname  = formatC(uthrname,digits=1,format="f")
    groups<-prune(gnmDend,thrd)
    write.table(groups,
                file=paste(paste("$outFile",thrname,sep="-"),"groups",sep="."),
                row.names=T,col.names=F,quote=F,append=F)
}
APPENDANI
    }
    elsif( $cutGroups ne "none" ) {
        print {$CLUSTANI} << "APPENDANI1";
thrd     = 100 - $cutGroups
thrname  = formatC($cutGroups,digits=1,format="f")
groups<-prune(gnmDend,thrd)
write.table(groups,
            file=paste(paste("$outFile",thrname,sep="-"),"groups",sep="."),
            row.names=T,col.names=F,quote=F,append=F)
APPENDANI1
    }
    close($CLUSTANI);
}

##### a subroutine to figure out file compression, if any, and
##### the proper program to open it
sub figureCompression {
    my $rootName = $_[0];
    $rootName =~ s{\.(gz|bz2|Z)$}{};
    my $fullName
        = ( -f "$rootName.gz" )  ? "$rootName.gz"
        : ( -f "$rootName.Z" )   ? "$rootName.Z"
        : ( -f "$rootName.bz2" ) ? "$rootName.bz2"
        : ( -f "$rootName" )     ? "$rootName"
        : "none";
    my $catProg
        = $fullName =~ m{\.(gz|Z)$} ? "gzip -qdc"
        : $fullName =~ m{\.bz2$}    ? "bzip2 -qdc"
        : "cat";
    if( $fullName eq "none" ) {
        return();
    }
    else {
        return("$catProg","$fullName");
    }
}

sub runR {
    my $rFile = $_[0];
    my $Rcmd = qq(R -f $rFile.R >& $rFile.log);
    print "      running R\n";
    chdir("$tmpDir");
    my $runR = qx($Rcmd);
    chdir("$cwd");
}

#### Cluster by Mash's Jaccard distance
sub clustMash {
    my($clustMethod,$subFile,$outFile,$cutGroups) = @_;
    my $rFile = "$tmpDir/$outFile.R";
    print "      producing R file:\n"
        . "$rFile\n";
    open( my $CLUSTMASH,">","$rFile" );
    print {$CLUSTMASH} << "CLUSTMASH";
$Rhead

############## import MASH data
mashtbl<-read.table("$subFile",header=T,sep="\\t")

##### obtain a vector of genome names to be able to name the
##### items in redundant group later
number.genome<-length(levels(as.factor(mashtbl[,2])))
genome.names<-mashtbl[1:number.genome,2]

##### define the distance matrix of using mashtbl
##### MASH (MASH = 3)
dx<-xpnd(mashtbl[,3])
rownames(dx)<-genome.names
colnames(dx)<-genome.names

##### the output file contains the distance matrix
xmat = as.matrix(dx)
write.table(xmat,
        bzfile("$root.dist.matrix.bz2",compression = 9),
        row.names=T,col.names=T,quote=F)

#### the output file contains the distance matrix as a table
xmat[upper.tri(xmat)] <- NA
x = melt(xmat, na.rm=T)
write.table(x, bzfile("$root.dist.tbl.bz2",compression = 9),
               quote=F,row.names=F,col.names=F,sep="\\t")

##### create a dendrogram using the agglomerative (bottom up, agnes)
##### or divisive (top-down, diana) and the dx dissimilarity matrix
dx<-as.dist(dx)
gnmCluster<-$procedure
gnmDend<-as.hclust(gnmCluster)
saveRDS(gnmDend,"${outFile}-cluster.Rds")

##### Output the tree in the Newick format
phy <- as.phylo(as.hclust(gnmCluster))
write.tree(phy, file="$outFile.tree")

##### calculate the divisive/agglomerative coefficients
coeff.gnm<-coef(gnmCluster)
coeff.gnm
coeff.print<-formatC(coeff.gnm,digits=3,format="f")
write(paste("$clustMethod", coeff.print),
      "$outFile.coeff",sep = "\\t", append=T)
CLUSTMASH
    if( $cutGroups eq "1:10" ) {
        print {$CLUSTMASH} << "APPENDMASH";
for ( i in 1:12 ){
    thrd<-(i + 4) * 0.005
    thrname = formatC(thrd,digits=3,format="f")
    groups<-prune(gnmDend,thrd)
    write.table(groups,
                file=paste(paste("$outFile",thrname,sep="-"),"groups",sep="."),
                row.names=T,col.names=F,quote=F,append=F)
}
APPENDMASH
    }
    elsif( $cutGroups ne "none" ) {
        print {$CLUSTMASH} << "APPENDMASH1";
thrd    = $cutGroups
thrname = formatC(thrd,digits=3,format="f")
groups<-prune(gnmDend,thrd)
write.table(groups,
            file=paste(paste("$outFile",thrname,sep="-"),"groups",sep="."),
            row.names=T,col.names=F,quote=F,append=F)
APPENDMASH1
    }
    close($CLUSTMASH);
}

sub fixMash {
    my $file2fix = $_[0];
    print "   cleaning:\n$file2fix\n";
    my ($open,$cfile) = figureCompression("$file2fix");
    my $cleanFl = $file2fix;
    $cleanFl =~ s{^\S+/}{};
    $cleanFl =~ s{\.\S+}{};
    $cleanFl = "$tmpDir/$cleanFl.clean.bz2";
    my %dist  = ();
    my %count = ();
    open( my $INMASH,"-|","$open $cfile" );
    while(<$INMASH>) {
        my($gnm1,$gnm2,$dist) = split;
        $gnm1 =~ s{^\S+/}{};
        $gnm2 =~ s{^\S+/}{};
        $gnm1 =~ s{\.\S+}{};
        $gnm2 =~ s{\.\S+}{};
        $count{"$gnm1"}++;
        $count{"$gnm2"}++;
        my $pair = join("\t",sort($gnm1,$gnm2));
        $dist{"$pair"} = $dist;
    }
    close($INMASH);
    my @genomes = sort keys %count;
    open( my $OUTMASH,"|-","bzip2 --best >$cleanFl" );
    print {$OUTMASH} join("\t","Genome1","Genome2","Distance"),"\n";
    while( my $gnm1 = shift @genomes ) {
        print {$OUTMASH} join("\t",$gnm1,$gnm1,"0"),"\n";
        for my $gnm2 ( @genomes ) {
            my $pair = join("\t",sort($gnm1,$gnm2));
            my $dist = exists $dist{"$pair"} ? $dist{"$pair"} : 1;
            print {$OUTMASH} join("\t",$gnm1,$gnm2,$dist),"\n";
        }
    }
    close($OUTMASH);
    return("$cleanFl");
}

sub fixDashing {
    my $file2fix = $_[0];
    my $fixN = $basis eq "dashing" ? "0" : "1";
    print "   cleaning (fix = $fixN):\n$file2fix\n";
    my ($open,$cfile) = figureCompression("$file2fix");
    my $cleanFl = $file2fix;
    $cleanFl =~ s{^\S+/}{};
    $cleanFl =~ s{\.\S+}{};
    $cleanFl = "$tmpDir/$cleanFl.clean.bz2";
    open( my $OUTDASH,"|-","bzip2 --best >$cleanFl" );
    open( my $INDASH,"-|","$open $cfile" );
    while(<$INDASH>) {
        if( m{^#} ) {
            s{^#\s*Names\s*}{}i;
            chomp;
            my @heading = ();
            for my $item ( split(/\t/,$_) ) {
                $item =~ s{\S+/}{};
                $item =~ s{\.\S+}{};
                push(@heading,$item);
            }
            print {$OUTDASH} join("\t",@heading),"\n";
        }
        else {
            my($gnm1,@items) = split;
            $gnm1 =~ s{^\S+/}{};
            $gnm1 =~ s{\.\S+}{};
            my @clean = ();
            for my $item ( @items ) {
                my $add = $item eq '' ? $fixN : $item;
                push(@clean,$add);
            }
            print {$OUTDASH} join("\t",$gnm1,@clean),"\n";
        }
    }
    close($INDASH);
    close($OUTDASH);
    return("$cleanFl");
}

#### Cluster Dashing's matrix (either mash distances or Jaccard indexes)
sub clustDashing {
    my($clustMethod,$subFile,$outFile,$cutGroups) = @_;
    my $fixMatrix
        = $basis eq "dashing"
        ? qq(diag(dashingtbl) <- 1\ndashingtbl = 1 - dashingtbl\n)
        : '';
    my $rFile = "$tmpDir/$outFile.R";
    print "      producing R file:\n"
        . "$rFile\n";
    open( my $CLUSTDASHING,">","$rFile" );
    print {$CLUSTDASHING} << "CLUSTDASHING";
$Rhead

############## import DASHING data
dashingtbl<-read.table("$subFile",sep="\\t")
$fixMatrix

##### the output file contains the distance matrix
xmat = as.matrix(dashingtbl)
write.table(xmat,
        bzfile("$root.dist.matrix.bz2",compression = 9),
        row.names=T,col.names=T,quote=F)

#### the output file contains the distance matrix as a table
xmat[upper.tri(xmat)] <- NA
x = melt(xmat, na.rm=T)
write.table(x, bzfile("$root.dist.tbl.bz2",compression = 9),
               quote=F,row.names=F,col.names=F,sep="\\t")

##### create a dendrogram using the agglomerative (bottom up, agnes)
##### or divisive (top-down, diana) and the dx dissimilarity matrix
dx<-as.dist(dashingtbl)
gnmCluster<-$procedure
gnmDend<-as.hclust(gnmCluster)
saveRDS(gnmDend,"${outFile}-cluster.Rds")

##### Output the tree in the Newick format
phy <- as.phylo(as.hclust(gnmCluster))
write.tree(phy, file="$outFile.tree")

##### calculate the divisive/agglomerative coefficients
coeff.gnm<-coef(gnmCluster)
coeff.gnm
coeff.print<-formatC(coeff.gnm,digits=3,format="f")
write(paste("$clustMethod", coeff.print),
      "$outFile.coeff",sep = "\\t", append=T)
CLUSTDASHING
    if( $cutGroups eq "1:10" ) {
        print {$CLUSTDASHING} << "APPENDDASHING";
for ( i in 1:15 ){
    thrd<-(i + 8) * 0.005
    thrname = formatC(thrd,digits=3,format="f")
    groups<-prune(gnmDend,thrd)
    write.table(groups,
                file=paste(paste("$outFile",thrname,sep="-"),"groups",sep="."),
                row.names=T,col.names=F,quote=F,append=F)
}
APPENDDASHING
    }
    elsif( $cutGroups ne "none" ) {
        print {$CLUSTDASHING} << "APPENDDASHING1";
thrd = $cutGroups
thrname = formatC(thrd,digits=3,format="f")
groups<-prune(gnmDend,thrd)
write.table(groups,
            file=paste(paste("$outFile",thrname,sep="-"),"groups",sep="."),
            row.names=T,col.names=F,quote=F,append=F)
APPENDDASHING1
    }
    close($CLUSTDASHING);
}

########################################################
########################################################
########################################################
#### Cluster domains by Jaccard distances
########################################################
########################################################
sub clustJaccard {
    my($clustMethod,$subFile,$outFile,$cutGroups,$weight) = @_;
    my $addLibrary = $weight eq "weighted" ? qq(library("philentropy")) : '';
    my $distance
        = $weight eq "weighted"
        ? 'distance(domaintbl,method="jaccard",use.row.names=TRUE)'
        : 'dist(domaintbl,method="binary")';
    my $rFile = "$tmpDir/$outFile.R";
    print "      producing R file:\n"
        . "$rFile\n";
    open( my $CLUSTJACCARD,">","$rFile" );
    print {$CLUSTJACCARD} << "CLUSTJACCARD";
$Rhead
$addLibrary
############## import the distance data
domaintbl<-read.table("$subFile", sep="\\t")

############## calculate the distance matrix
dx <- $distance

#### the output file contains the distance matrix
xmat = as.matrix(dx)
write.table(xmat,
        bzfile("${root}-$weight.dist.matrix.bz2",compression = 9),
        row.names=T,col.names=T,quote=F)

#### the output file contains the distance matrix as a table
xmat[upper.tri(xmat)] <- NA
x = melt(xmat, na.rm=T)
write.table(x, bzfile("${root}-$weight.dist.tbl.bz2",compression = 9),
               quote=F,row.names=F,col.names=F,sep="\\t")

## create a dendrogram using the agglomerative (bottom up, agnes)
## or divisive (top-down, diana) and the dx dissimilarity matrix
gnmCluster<-$procedure
gnmDend<-as.hclust(gnmCluster)
saveRDS(gnmDend,"${outFile}-cluster.Rds")

## Output the tree in the Newick format
phy <- as.phylo(as.hclust(gnmCluster))
write.tree(phy, file="$outFile.tree")

## calculate the divisive/agglomerative coefficients
coeff.gnm<-coef(gnmCluster)
coeff.gnm
coeff.print<-formatC(coeff.gnm,digits=3,format="f")
write(paste("$clustMethod", coeff.print),
      "$outFile.coeff",sep = "\\t", append=T)
CLUSTJACCARD
    if( $cutGroups eq "1:10" ) {
        print {$CLUSTJACCARD} << "APPENDJB";
for ( i in 1:15 ){
    thrd<-i * 0.01
    thrname = formatC(thrd,digits=3,format="f")
    groups<-prune(gnmDend,thrd)
    write.table(groups,
                file=paste(paste("$outFile",thrname,sep="-"),"groups",sep="."),
                row.names=T,col.names=F,quote=F,append=F)
}
APPENDJB
    }
    elsif( $cutGroups ne "none" ) {
        print {$CLUSTJACCARD} << "APPENDJB1";
thrd = $cutGroups
thrname = formatC(thrd,digits=3,format="f")
groups<-prune(gnmDend,thrd)
write.table(groups,
            file=paste(paste("$outFile",thrname,sep="-"),"groups",sep="."),
            row.names=T,col.names=F,quote=F,append=F)
APPENDJB1
    }
    close($CLUSTJACCARD);
}
