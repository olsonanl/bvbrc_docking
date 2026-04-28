#
# App wrapper for the diffdock-based docking pipeline.
#
# Docking pipeline requires an invocation of run_local_docking for each PDB file and set of ligands.
#
# The output folder includes a copy of the PDB file used to run the docking, as well as a log
# of the underlying diffdock invocation. Each ligand has an output folder labeled by
# the tag given to the ligand. For the input option where a user-defined list of smiles strings
# is provided, the app will generate ids of the form ligand-XX.
#
# We also generate a report with a tabular summary of the entire run, including links to the
# generated ligand-in-structure PDBs.
#

use Carp::Always;
use Bio::KBase::AppService::AppScript;
use File::Slurp;
use IPC::Run;
use Cwd qw(abs_path getcwd);
use File::Path qw(rmtree make_path);
use strict;
use Data::Dumper;
use File::Basename;
use File::Temp;
use JSON::XS;
use Getopt::Long::Descriptive;
use DockingCompute;

my $docking_compute = new DockingCompute();

my $app = Bio::KBase::AppService::AppScript->new(sub { $docking_compute->run(@_); },
						 sub { $docking_compute->preflight(@_); });

$app->run(\@ARGV);

