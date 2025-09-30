#! /usr/bin/env perl
#
# The StabilityPrediction application.
#

use Bio::KBase::AppService::AppScript;
use strict;
use Data::Dumper;
use Cwd;
use File::Basename;
use File::Path qw(remove_tree);

my $script = Bio::KBase::AppService::AppScript->new(\&stabPred, \&preflight);

$script->run(\@ARGV);

# Count standard AA residues (unique chain:resseq:icode) from ATOM records.
# - Ignores HETATM, waters, ligands.
# - Keeps only altLoc = ' ' or 'A' (primary conformer).
# - Collapses duplicates across atoms in the same residue.
sub count_pdb_residues {
    my ($pdb_path) = @_;
    open my $fh, '<', $pdb_path or die "Cannot open $pdb_path: $!";

    my %std = map { $_ => 1 } qw(
        ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET
        PHE PRO SER THR TRP TYR VAL
    );

    my %seen_res;          # key: "chain:resseq:icode" -> 1
    my %first_alt_for_res; # track if we've already accepted a primary altLoc for a residue

    while (my $line = <$fh>) {
        next unless length($line) >= 27;           # ensure columns exist
        next unless substr($line, 0, 4) eq 'ATOM'; # only ATOM records

        my $altLoc  = substr($line, 16, 1);                  # col 17
        my $resName = substr($line, 17, 3); $resName =~ s/\s+//g;  # cols 18-20
        my $chainID = substr($line, 21, 1); $chainID =~ s/\s+//g;  # col 22
        my $resSeq  = substr($line, 22, 4); $resSeq  =~ s/\s+//g;  # cols 23-26
        my $iCode   = substr($line, 26, 1); $iCode   =~ s/\s+//g;  # col 27

        next unless $std{$resName};  # only standard amino acids
        my $key = join(':', $chainID, $resSeq, $iCode);

        # Accept only primary conformer per residue (blank or 'A')
        # If we've already accepted a primary for this residue, skip others.
        if ($altLoc ne ' ' && $altLoc ne 'A') {
            next unless $first_alt_for_res{$key};  # only allow if primary already seen (we ignore non-primary)
            next;                                  # in practice, we just skip non-primary entirely
        } else {
            $first_alt_for_res{$key} ||= 1;        # mark primary as seen
        }

        $seen_res{$key} = 1;
    }
    close $fh;

    return scalar keys %seen_res;
}

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;

    # set temp directory
    my $tmp_dir = $params->{_tmpdir};
    $tmp_dir //= getcwd . "/tmp.$$";
    mkdir($tmp_dir) or die "Failed to make temp directory\n";

    # initialize PDB file to use
    my $pdbFile;

    # If no PDB file and no PDB ID is provided, kill process as we 
    # need to compute on something.  
    if ( $params->{pdb_id} eq "" and $params->{pdb} eq "" ) {
        die "PDB ID or PDB File must be provided";
    }

    # if both PDB file and PDB ID are provided, kill process as we
    # compute on one or the other.
    # currently not planning on supporting batch PDB files
    if ( $params->{pdb_id} ne "" and $params->{pdb} ne "" ) {
        die "PDB ID and PDB File both provided, must provide one or the other";
    }

    # if a PDB ID is provided
    #   Download the PDB file from PDBBank
    #   Set the PDB file
    if ( $params->{pdb_id} ne "" ) {
        $pdbFile = $tmp_dir . "/" . $params->{pdb_id} . ".pdb";
        # initialize the return code to 0
        my $rc = 0;
        # download file
        $rc = system("curl", "-o", $pdbFile, "https://files.rcsb.org/download/" . $params->{pdb_id} . ".pdb");
        # if download failed, output error and exit
        die "Failure downloading PDB file" if $rc != 0;
    # otherwise, set the PDB file from the one provided from 
    #   the workspace
    } else {
        $pdbFile = $tmp_dir . "/" . basename($params -> {pdb});
        my $rc = system("p3-cp", "ws:" . $params->{pdb}, $pdbFile);
        die "Failure copying PDB file" if $rc != 0;
    }

    # Preflight requirements based on double the memory and runtime
    #   of various tests
    # Epistatic runs take a LOT longer to run than single or additive
    # For epistatic runs, set runtime to 1 day
    if($params->{mode} eq "epistatic") {
        # base runtime and memory on the number of residues in the PDB
        # added a few buffers for better fitment
        my $res = count_pdb_residues($pdbFile);

        # runtime and memory requirements based on timed test runs
        my $runtime = int((0.0248 * $res ** 2 - 10.294 * $res + 4956.8) * 1.2);
        my $memReq = int((152020 * $res) * 1.2 / 1000000 * 1.2);

        # set CPU, memory requirements, and runtime
        # CPU = 8 is optimal based on rough testing 
        #   for epistatic mode
        return { cpu => 8, memory => $memReq . "G", runtime => $runtime };
    }
    elsif($params->{mode} eq "additive") {
        # base runtime and memory on the number of residues in the PDB
        my $res = count_pdb_residues($pdbFile);

        # runtime and memory requirements based on test runs
        my $runtime = int((3.9802 * $res) * 1.2);
        my $memReq = int((137745 * $res) * 1.2);

        # set CPU, memory requirements, and runtime
        return { cpu => 1, memory => $memReq . "G", runtime => $runtime };
    }
    else {
        my $res = count_pdb_residues($pdbFile);
        my $memReq = int((639.77 * $res + 567236) * 1.2);
        # set CPU, memory requirements, and runtime
        # Runtime typically less than 1 minute, but give it a buffer
        #   that is still reasonably small
        # CPU = 1 to minimize CPU requirements
        return { cpu => 1, memory => $memReq . "G", runtime => 600 };
    }
}

sub stabPred
{
    my($app, $app_def, $raw_params, $params) = @_;

    # set temp directory
    my $tmp_dir = $params->{_tmpdir};
    $tmp_dir //= getcwd . "/tmp.$$";
    mkdir($tmp_dir) or die "Failed to make temp directory\n";

    # initialize PDB file to use
    my $pdbFile;

    # If no PDB file and no PDB ID is provided, kill process as we 
    # need to compute on something.  
    if ( $params->{pdb_id} eq "" and $params->{pdb} eq "" ) {
        die "PDB ID or PDB File must be provided";
    }

    # if both PDB file and PDB ID are provided, kill process as we
    # compute on one or the other.
    # currently not planning on supporting batch PDB files
    if ( $params->{pdb_id} ne "" and $params->{pdb} ne "" ) {
        die "PDB ID and PDB File both provided, must provide one or the other";
    }

    # if a PDB ID is provided
    #   Set the PDB file name from downloaded file in preflight
    if ( $params->{pdb_id} ne "" ) {
        $pdbFile = $tmp_dir . "/" . $params->{pdb_id} . ".pdb";
        my $rc = 0;
        # download file
        $rc = system("curl", "-o", $pdbFile, "https://files.rcsb.org/download/" . $params->{pdb_id} . ".pdb");
        # if download failed, output error and exit
        die "Failure downloading PDB file" if $rc != 0;
    # otherwise, set the PDB file from the one provided from 
    #   the workspace
    } else {
        $pdbFile = $tmp_dir . "/" . basename($params -> {pdb});
        my $rc = system("p3-cp", "ws:" . $params->{pdb}, $pdbFile);
        die "Failure copying PDB file" if $rc != 0;
    }

    # define number of threads
    my $threads = $ENV{P3_ALLOCATED_CPU} // 1;

    # initialize the return code to 0
    my $rc = 0;
    my $cmnd = "";
    # if chains is specified, then run the command with the 
    # chains options
    if ( $params->{chains} ne "" ) {
        # if ss_penalty is true, then add ss_penalty option
        if ( $params->{ss_penalty} eq 1 ) {
            $rc = system("ThermoMPNN-D", "--mode", $params->{mode}, "--pdb", $pdbFile, "--batch_size", $params->{batch_size}, "--out", "thermompnnd", "--chains", $params->{chains}, "--threshold", $params->{threshold}, " --distance ", $params->{distance}, " --ss_penalty ", "--threads", $threads);
            $cmnd = "ThermoMPNN-D" . " --mode " . $params->{mode} . " --pdb " . $pdbFile . " --batch_size " . $params->{batch_size} . " --out " . "thermompnnd" . " --chains " . $params->{chains} . " --threshold " . $params->{threshold} . " --distance " . $params->{distance} . " --ss_penalty " . " --threads " . $threads;
        } else {
            $rc = system("ThermoMPNN-D", "--mode", $params->{mode}, "--pdb", $pdbFile, "--batch_size", $params->{batch_size}, "--out", "thermompnnd", "--chains", $params->{chains}, "--threshold", $params->{threshold}, "--distance", $params->{distance}, "--threads", $threads);
            $cmnd = "ThermoMPNN-D" . " --mode " . $params->{mode} . " --pdb " . $pdbFile . " --batch_size " . $params->{batch_size} . " --out " . "thermompnnd" . " --chains " . $params->{chains} . " --threshold " . $params->{threshold} . " --distance " . $params->{distance} . " --ss_penalty " . " --threads " . $threads;
        }
    }
    # otherwise run the command without the chains option
    # running with "None" or "" causes it to crash as it defaults to
    # looking for A or some other chain.
    else {
        # if ss_penalty is true, then add ss_penalty option
        if ( $params->{ss_penalty} eq 1 ) {
            $rc = system("ThermoMPNN-D", "--mode", $params->{mode}, "--pdb", $pdbFile, "--batch_size", $params->{batch_size}, "--out", "thermompnnd", "--threshold", $params->{threshold}, "--distance", $params->{distance}, "--ss_penalty", "--threads", $threads);
            $cmnd = "ThermoMPNN-D" . " --mode " . $params->{mode} . " --pdb " . $pdbFile . " --batch_size " . $params->{batch_size} . " --out " . " thermompnnd " . " --threshold " . $params->{threshold} . " --distance " . $params->{distance} . " --ss_penalty " . " --threads " . $threads;
        } else {
            $rc = system("ThermoMPNN-D", "--mode", $params->{mode}, "--pdb", $pdbFile, "--batch_size", $params->{batch_size}, "--out", "thermompnnd", "--threshold", $params->{threshold}, "--distance", $params->{distance}, "--threads", $threads);
            $cmnd = "ThermoMPNN-D" . " --mode " . $params->{mode} . " --pdb " . $pdbFile . " --batch_size " . $params->{batch_size} . " --out " . "thermompnnd" . " --threshold " . $params->{threshold} . " --distance " . $params->{distance} . " --threads " . $threads;
        }
    }

    # if the process fails, output error message and die
    die "Failure running ThermoMPNN-D: " . $cmnd  if $rc != 0;

    # get the result folder from the app specs
    my $folder = $app->result_folder();
    # output the thermompnnd.csv file output from the process to
    # the results folder, set the file type to csv
    # save file to file arguments
    #   file save file from (process output)
    #   file metadata to store
    #   file to save output to (workspace)
    #   file type
    #   use shock to store the file
    $app->workspace->save_file_to_file("thermompnnd.csv", {}, "$folder/thermompnnd.csv", "csv", 1);
    $rc = remove_tree($tmp_dir);
}
