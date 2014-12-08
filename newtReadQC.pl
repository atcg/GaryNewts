#!/usr/bin/perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Cwd;
use Parallel::ForkManager;

my $help = 0;
my $readsDir;
my $adaptersDir;
my $outDir;
my $logFile;
my $threadsMax = 4;
my $trim; # If set, will run Trimmomatic
my $filter; # If set, it will remove reads that were filtered by CASAVA
my $join; # If set, will run fastq-join to merge overlapping paired end reads
my $cat; # If set, will combine the singleton and joined reads into a single file


GetOptions  ("reads=s"         => \$readsDir,
             "adapters=s"      => \$adaptersDir,
             "out=s"           => \$outDir,
             "log=s"           => \$logFile,
             "threads=i"       => \$threadsMax,
             "filter"          => \$filter,
             "trim"            => \$trim,
             "join"            => \$join,
             "cat"             => \$cat,
             "help|man" => \$help) || pod2usage(2);

if (!$readsDir or !$adaptersDir or !$outDir or $help) {
    pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1);
}

my $originalDir = getcwd();

my @samples = ("Gary1", "Gary2", "HBS124214", "HBS124405", "HBS124407");
open(my $logFH, ">", $logFile) or die "Couldn't open log file $logFile for writing: $!\n";


# Make sure we have the right directories for output files
my $startingDir = getcwd();
unless (-d $outDir) {
    mkdir $outDir or die "Couldn't make output directory $outDir: $!\n";
}

my $fqjDir = $outDir . "/fastq-join";
unless (-d $fqjDir) {
    mkdir $fqjDir;
}

my $trimmomaticDir = $outDir . "/trimmomatic";
unless (-d $trimmomaticDir) {
    mkdir $trimmomaticDir or die "Couldn't make directory $trimmomaticDir: $!\n";
}



########## Remove the :Y: files that were designated as bad by CASAVA ##########
if ($filter) {
    print $logFH "--------------------------------------------------\n";
    print $logFH "Removing reads that were filtered by CASAVA\n";
    print $logFH "--------------------------------------------------\n";
    my $filterForkManager = Parallel::ForkManager->new($threadsMax);
    foreach my $tort (@samples) {
        $filterForkManager->start and next;
        my $R1File = $readsDir . $tort . "_R1.fastq.gz";
        my $R2File = $readsDir . $tort . "_R2.fastq.gz";
        my $R1outFile = $readsDir . $tort . "_R1_Ns.fastq.gz";
        my $R2outFile = $readsDir . $tort . "_R2_Ns.fastq.gz";
        system("gunzip -c $R1File | fastq_illumina_filter -vN | gzip > $R1outFile");
        system("gunzip -c $R2File | fastq_illumina_filter -vN | gzip > $R2outFile");
        $filterForkManager->finish;
    }
    $filterForkManager->wait_all_children;
    print $logFH "--------------------------------------------------\n";
    print $logFH "Finished removing reads that were filtered by CASAVA\n\n";
    print $logFH "--------------------------------------------------\n";
}


########## Generate Trimmomatic commands and run Trimmomatic ##########
print $logFH "--------------------------------------------------\n";
print $logFH "Generating Trimmomatic commands for all read files\n";
print $logFH "--------------------------------------------------\n";
my @trimmomaticCommands;
foreach my $tort (@samples) {
    my $R1File = $readsDir . "$tort" . "_R1_Ns.fastq.gz";
    my $R2File = $readsDir . "$tort" . "_R2_Ns.fastq.gz";
    my $adaptersFile = $adaptersDir . $tort . "_adapters.fasta";
    my $R1OutFilePaired = $trimmomaticDir . "/$tort" . "_R1p_Ns_trim.fastq.gz";
    my $R2OutFilePaired = $trimmomaticDir . "/$tort" . "_R2p_Ns_trim.fastq.gz";
    my $R1OutFileSingles = $trimmomaticDir . "/$tort" . "_R1s_Ns_trim.fastq.gz";
    my $R2OutFileSingles = $trimmomaticDir . "/$tort" . "_R2s_Ns_trim.fastq.gz";

    push (@trimmomaticCommands, "java -Xmx8g -jar /home/evan/bin/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 2 -phred33 $R1File $R2File $R1OutFilePaired $R1OutFileSingles $R2OutFilePaired $R2OutFileSingles ILLUMINACLIP:$adaptersFile:2:30:10 LEADING:5 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:40");    
}
print $logFH "--------------------------------------------------\n";
print $logFH "Finished generating trimmomatic commands for all tortoises\n";
print $logFH "--------------------------------------------------\n";


# Now run Trimmomatic
if ($trim) {
    print $logFH "--------------------------------------------------\n";
    print $logFH "Running all trimmomatic commands\n";
    print $logFH "--------------------------------------------------\n";
    my $counter = 0;
    my $trimThreads = $threadsMax / 2; # We're using two threads for every Trimmomatic process
    my $trimForkManager = new Parallel::ForkManager($threadsMax);
    foreach my $trimCommand (@trimmomaticCommands) {
        $counter++;
        print $logFH "--------------------------------------------------\n";
        print $logFH "Trimmomatic command $counter: \n\t"; # Indent the next line to make it easier to find the commands in the text
        print $logFH $trimCommand . "\n";
        sleep 3;
        print "\n";
        $trimForkManager->start and next;
        print "\n";
        system("$trimCommand");
        print "Finished running the following:\n\t$trimCommand\n\n";
        $trimForkManager->finish;
    }
    $trimForkManager->wait_all_children;
    print $logFH "--------------------------------------------------\n";
    print $logFH "Finished running all trimmomatic commands\n";
    print $logFH "--------------------------------------------------\n\n";
}


# Now run fastq-join to merge overlapping reads
if ($join) {
    print $logFH "--------------------------------------------------\n";
    print $logFH "Running fastq-join on all of the samples\n";
    print $logFH "--------------------------------------------------\n";
    my $fqjThreadsMax = $threadsMax / 4; # When using gzipped files each process uses about 4 threads
    my $joinForkManager = Parallel::ForkManager->new($fqjThreadsMax);
    foreach my $tort (@samples) {
        $joinForkManager->start and next;
        my $R1File = $trimmomaticDir . "/$tort" . "_R1p_Ns_trim.fastq.gz";
        my $R2File = $trimmomaticDir . "/$tort" . "_R2p_Ns_trim.fastq.gz";
        my $outPrefix = $fqjDir . "/$tort" . "_Ns_trim_fqj.%.fastq.gz";
        system("fastq-join -v ' ' $R1File $R2File -o $outPrefix");
        print $logFH "Finished running fastq-join on $tort\n";
        $joinForkManager->finish;
    }
    $joinForkManager->wait_all_children;
    print $logFH "--------------------------------------------------\n";
    print $logFH "Finished running fastq-join on all samples\n";
    print $logFH "--------------------------------------------------\n\n";
}


# Now combine all the singleton and joined reads into a single-end file
if ($cat) {
    print $logFH "--------------------------------------------------\n";
    print $logFH "Combining singleton and joined reads into single-end files\n";
    print $logFH "--------------------------------------------------\n";
    my $catForkManager = Parallel::ForkManager->new($threadsMax);
    foreach my $tort (@samples) {
        my $singleEndsFile = "$fqjDir/$tort" . "_Ns_trim_fqj_joinedAndSingles.fastq.gz";
        my $catCommand = "cat " . "$fqjDir/$tort" . "_Ns_trim_fqj.join.fastq.gz " . "$trimmomaticDir/$tort" . "_R1s_Ns_trim.fastq.gz " . "$trimmomaticDir/$tort" . "_R2s_Ns_trim.fastq.gz > $singleEndsFile";
        print "Running the following concatenation command:\n$catCommand\n";
        $catForkManager->start and next;
        system($catCommand);
        $catForkManager->finish;
    }
    $catForkManager->wait_all_children;
    print $logFH "--------------------------------------------------\n";
    print $logFH "Finished combining singleton and joined reads files\n";
    print $logFH "--------------------------------------------------\n";
}











#Documentation
__END__

=head1 NAME

newtReadQC.pl

=head1 SYNOPSIS 

perl newtReadQC.pl --reads <file> --adapters <file> --out <outputDirectory> --log <logfile.txt> --trim --threads 6 --reference <path/to/bwa/reference/genome/index.fasta>

 Options:
   -reads=s           Directory with raw reads in gzipped fastq format
   -adapters=s        Directory with adapters fasta files
   -out=s             Name of output directory (it will be created)
   -log=s             Name of logfile to print output (you will probably also want
                      to capture STDERR manually)
   -filter
   -trim              (flag) Perform read trimming using Trimmomatic and join with fastq-join
   -join              Run fastq-join
   -cat               Combine the singleton and joined read files into a single file
   -threads=i         Threads to use for multithreading (default=4)
   -help|man          Prints out documentation


=head1 DESCRIPTION

This was written for the purposes of QCing HiSeq reads for a desert tortoise
project

=cut