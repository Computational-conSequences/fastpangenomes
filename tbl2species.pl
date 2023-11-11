#!/usr/bin/env perl

use Getopt::Long;
use strict;

my $speciesTbl  = '';
my $taxonomyTbl = '';
my @distTbls    = ();
my $resultsDir  = 'SameTaxaTables';
my $defMix      = 'T';
my $mixIt       = $defMix;

my $ownName = $0;
$ownName =~ s{.*/}{};
my $helpBrief
    = qq(about:\n)
    . qq(    This program adds a "Same" column with species, genera, or\n)
    . qq(    family labels for genome identifier pairs based on NCBI's\n)
    . qq(    taxonomy\n\n)
    . qq(usage:\n)
    . qq(    $ownName -s <species-table> -t <taxonomy-table> -c <cluster-files(s)> \n\n)
    . qq(examples:\n)
    . qq(    $ownName -s DataSets/Enterobacteriaceae.species -t Complete.taxonomy.tsv.bz2 -c Clusters/*.dist.tbl.bz2\n\n)
    . qq(options:\n)
    . qq(   -s table with species to genome map [Complete.groups.bz2],\n)
    . qq(       required\n)
    . qq(   -t table with full genome taxonomy [Complete.taxonomy.bz2],\n)
    . qq(       required\n)
    . qq(   -c pairwise distance tables [ANI.dist.tbl.bz2], obtained with\n)
    . qq(       clusterGenomes.pl, required\n)
    . qq(   -m mix distances into a single table [T|F], default $defMix\n)
    . qq(   -o output directory, default '$resultsDir'\n)
    . qq(\n)
    ;

my $options = GetOptions(
    "s=s"    => \$speciesTbl,
    "t=s"    => \$taxonomyTbl,
    "c=s{,}" => \@distTbls,
    "o=s"    => \$resultsDir,
    "m=s"    => \$mixIt,
) or die $helpBrief;

my $outputFl = $speciesTbl;
$outputFl =~ s{\S+/}{};
$outputFl =~ s{\.\S+}{};

if( !$speciesTbl && $taxonomyTbl && !@distTbls ) {
    die $helpBrief;
}
elsif( !$speciesTbl ) {
    die "I need a genome ID to species table (-s)\n$helpBrief";
}
elsif( !$taxonomyTbl ) {
    die "I need a full genome taxonomy table (-t)\n$helpBrief";
}
elsif( !@distTbls ) {
    die "I need pairwise distance tables (-c)\n$helpBrief";
}

my @trueTbls = sort checkFl(uniq(@distTbls));
my $cFiles = @trueTbls;
if( $cFiles > 0 ) {
    print "working with $cFiles files\n";
}
else {
    die "the files:" . join("\n",uniq(@distTbls)) ."\ndo not exist\n";
}
$mixIt
    = $cFiles == 1          ? 'F'
    : $mixIt =~ m{^(T|F)$}i ? uc($1)
    : $defMix;

print "learning species and genera from $speciesTbl\n";
my( $rSpecies,$rGenus,$rFamily ) = learnTaxa("$speciesTbl","$taxonomyTbl");
my $countSpecies = keys %{$rSpecies};
my $countGenus   = keys %{$rGenus};
my $countFamily  = keys %{$rFamily};
print "learned species/genera/families for $countSpecies/$countGenus/$countFamily genomes\n";

unless( -d "$resultsDir" ) {
    mkdir("$resultsDir");
}

if( $mixIt eq 'T' ) {
    print "mixing distances and adding columns to distance tables\n";
    my @heads = ();
    my @dists = ();
    push(@heads,naked($trueTbls[0]));
    my($rDist,$rClass,$rName) = addTaxa(shift @trueTbls,'Y');
    for my $distTbl ( @trueTbls ) {
        push(@heads,naked($distTbl));
        addTaxa($distTbl,'N',$rDist);
    }
    my $outFile = "$resultsDir/$outputFl.dist.bz2";
    my $tmpFile = $outFile . ".tmp";
    print "printing $outFile\n";
    my $cDists   = @heads;
    my $problems = 0;
    my $printed  = 0;
    open( my $OUTFL,"|-","bzip2 --best > $tmpFile" );
    print {$OUTFL} join("\t","Genome1","Genome2",@heads,"Same","Name"),"\n";
    for my $pair ( sort keys %{ $rDist } ) {
        my $cCheck = split(/\t/,$rDist->{"$pair"});
        if( $cCheck eq $cDists ) {
            print {$OUTFL} join("\t",
                                $pair,
                                $rDist->{"$pair"},
                                $rClass->{"$pair"},
                                $rName->{"$pair"}
                            ),"\n";
            $printed++;
        }
        else {
            $problems++;
        }
    }
    close($OUTFL);
    if( $printed > 0 ) {
        print "printed $printed pairs with $problems problematic pairs\n";
        rename($tmpFile,$outFile);
    }
    else {
        print "no pairs printed\n";
        unlink($tmpFile);
    }
}
else {
    print "adding columns to distance tables\n";
    for my $distTbl ( @trueTbls ) {
        print "   working with $distTbl\n";
        my $outFile = $distTbl;
        $outFile =~ s{\S+/}{};
        $outFile = $resultsDir . "/" . $outFile;
        my( $printed,$rejected )
            = addTaxa($distTbl,$outFile,$rSpecies,$rGenus);
        print "      added columns for $printed genome pairs\n"
            . "      rejected $rejected pairs for lack of a species name\n";
    }
}

print "        $ownName done\n\n";

#################################################################
#################################################################
####################### sub routines ############################
#################################################################
#################################################################

sub learnTaxa {
    my ($spFl,$txFl) = @_;
    my %species = ();
    my %genus   = ();
    my $open
        = $spFl =~ m{\.bz2$}    ? "bzip2 -qdc"
        : $spFl =~ m{\.(gz|Z)$} ? "gzip -qdc"
        : "cat";
    open( my $SPFL,"-|","$open $spFl" ) or die "problems with $spFl $!";
  READSPFL:
    while(<$SPFL>) {
        next READSPFL if m{^(Species|#)};
        chomp;
        my( $species,$count,@gcfs ) = split(/\t/,$_);
        my $genus = $species;
        $genus =~ s{^Candidatus\s+}{};
        $genus =~ s{\s+[a-z]+$}{};
        $genus =~ s{\s+[a-z]+$}{};
        for my $gcf ( @gcfs ) {
            $species{"$gcf"} = $species;
            $genus{"$gcf"}   = $genus;
        }
    }
    close($SPFL);
    my %family = ();
    my $opent
        = $txFl =~ m{\.bz2$}    ? "bzip2 -qdc"
        : $txFl =~ m{\.(gz|Z)$} ? "gzip -qdc"
        : "cat";
    open( my $TXFL,"-|","$opent $txFl" ) or die "problems with $txFl $!";
  READTX:
    while(<$TXFL>) {
        next READTX if m{^GenomeID};
        chomp;
        my($genomeid,$superkingdom,$phylum,$class,
           $order,$family,$genus,$strain) = split(/\t/,$_);
        if( $family ne "NA" && exists $species{"$genomeid"}) {
            $family{"$genomeid"} = $family;
        }
    }
    close($TXFL);
    return(\%species,\%genus,\%family);
}

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}

sub checkFl {
    my @test = @_;
    my $initial = @test;
    my @found   = ();
    my @missing = ();
    for my $distTbl ( @test ) {
        if( -f $distTbl ) {
            push(@found,$distTbl);
        }
        else {
            push(@missing,$distTbl);
        }
    }
    my $missing = @missing;
    my $found   = @found;
    if( $missing > 0 ) {
        print "$missing missing files:\n",
            join("\n",@missing),"\n";
    }
    if( $found > 0 ) {
        print "found $found files\n";
        return(@found);
    }
    else {
        return();
    }
}

sub addTaxa {
    my($table,$init,$rDist) = @_;
    print "extracting data from $table\n";
    my $printed   = 0;
    my $nospecies = 0;
    my $rejected  = 0;
    my %same = () if( $init eq 'Y' );
    my %name = () if( $init eq 'Y' );
    my %dist = ();
    ######### begin if not mixing ##########
    my $outFile = '';
    my $tmpFile = '';
    my $OUTFL;
    if( $mixIt eq 'F' ) {
        $outFile = $resultsDir . "/" . naked($table) . ".dist.bz2";
        $tmpFile = $outFile . ".tmp";
        open( $OUTFL,"|-","bzip2 --best > $tmpFile" );
        print {$OUTFL}
            join("\t","Genome1","Genome2",naked($table),"Same","Name"),"\n";
    }
    ######### end if not mixing   ##########
    open( my $DISTFL,"-|","bzip2 -qdc $table" );
  READDTBLC:
    while(<$DISTFL>) {
        my( $gcf1,$gcf2,$dist ) = split;
        next READDTBLC if( $gcf1 eq $gcf2 );
        if( !exists $rSpecies->{"$gcf1"}
            || !exists $rSpecies->{"$gcf2"} ) {
            $nospecies++;
            $rejected++;
            next READDTBLC;
        }
        if( !exists $rFamily->{"$gcf1"}
            || !exists $rFamily->{"$gcf2"} ) {
            $rejected++;
            next READDTBLC;
        }
        if( $rFamily->{"$gcf1"} ne $rFamily->{"$gcf2"} ) {
            $rejected++;
            next READDTBLC;
        }
        my $pair = join("\t",sort $gcf1,$gcf2);
        my($same,$name)
            = ( $rSpecies->{"$gcf1"} eq $rSpecies->{"$gcf2"} )
            ? ( "Species",$rSpecies->{"$gcf1"} )
            : ( $rGenus->{"$gcf1"} eq $rGenus->{"$gcf2"} )
            ? ( "Genus",$rGenus->{"$gcf1"} )
            : ( "Family",$rFamily->{"$gcf1"} );
        if( $mixIt eq 'T' ) {
            if( $init eq 'Y' ) {
                $dist{"$pair"} = $dist;
                $same{"$pair"} = $same;
                $name{"$pair"} = $name;
            }
            else {
                if( exists $rDist->{"$pair"} ) {
                    $rDist->{"$pair"} .= "\t" . $dist;
                }
            }
        }
        else { ### $mixIt eq 'F'
            print {$OUTFL}
                join("\t",$pair,$dist,$same,$name),"\n";
        }
        $printed++;
    }
    close($DISTFL);
    if( $mixIt eq 'F' ) {
        close($OUTFL);
        if( $printed > 0 ) {
            rename($tmpFile,$outFile);
        }
        else {
            print "    no printed pairs\n";
            unlink($tmpFile);
        }
        return($printed,$rejected);
    }
    elsif( $init eq 'Y' ) {
        print "   distances initialized\n";
        print "   added columns for $printed genome pairs\n"
            . "   rejected $rejected pairs for taxonomic problems\n"
            . "   of those $nospecies pairs lack a species name\n";
        return(\%dist,\%same,\%name);
    }
}

sub naked {
    my $clean = $_[0];
    $clean =~ s{^\S+/}{};
    $clean =~ s{\.\S+$}{};
    return($clean);
}
