package DockingCompute;

use File::Slurp;
use IPC::Run;
use Cwd qw(abs_path getcwd);
use File::Copy;
use File::Copy::Recursive qw(dircopy);
use File::Path qw(rmtree make_path);
use File::Spec;
use strict;
use Data::Dumper;
use File::Basename;
use File::Temp;
use JSON;
use Text::CSV qw(csv);
use Getopt::Long::Descriptive;
use P3DataAPI;
use YAML;
use Template;
use Module::Metadata;
use warnings;
use strict;

use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(api smiles_list debug work_dir staging_dir output_dir app params
			     ligand_name));

sub new
{
    my($class) = @_;
    my $self = {};
    bless $self, $class;

    $self->{api} = P3DataAPI->new();
    
    return $self;
}

sub run
{
    my($self, $app, $app_def, $raw_params, $params) = @_;

    #
    # Set up work and staging directories.
    # 
    $self->app($app);
    $self->params($params);

    my $tmp_dir = $params->{_tmpdir};
    $tmp_dir //= getcwd . "/tmp.$$";
    
    my $work_dir = "$tmp_dir/work";
    my $staging_dir = "$tmp_dir/staging";
    my $output_dir = "$tmp_dir/output";

    $self->{work_dir} = $work_dir;
    $self->{staging_dir} = $staging_dir;
    $self->{output_dir} = $output_dir;

    make_path($work_dir, $staging_dir, $output_dir);

    #
    # Find the lingands
    #
    my $ligand_file;
    if ($params->{ligand_library_type} eq 'named_library')
    {
        if ($params->{ligand_named_library} eq 'approved-drugs')
        {
	    $ligand_file = $self->load_ligand_library("/vol/bvbrc/production/application-backend/bvbrc_docking/drugbank_approved.txt");
        }
        elsif ($params->{ligand_named_library} eq 'experimental_drugs')
        {
	    $ligand_file = $self->load_ligand_library("/vol/bvbrc/production/application-backend/bvbrc_docking/drugbank_exp_inv.txt");
        }
        elsif ($params->{ligand_named_library} eq 'test')
        {
        $ligand_file = $self->load_ligand_library("/vol/bvbrc/production/application-backend/bvbrc_docking/test.txt");
        }
        elsif ($params->{ligand_named_library} eq 'small_db')
        {
        $ligand_file = $self->load_ligand_library("/vol/bvbrc/production/application-backend/bvbrc_docking/small_db.txt");
        }
        else
        {
        die "Unknown ligand library type selected $params->{ligand_library_type}";
        }

    }
    elsif ($params->{ligand_library_type} eq 'smiles_list')
    {
	$ligand_file = $self->load_ligand_smiles($params->{ligand_smiles_list});
    }
    elsif ($params->{ligand_library_type} eq 'ws_file')
    {
	$ligand_file = $self->load_ligand_ws_file($params->{ligand_ws_file});
    }

    #
    # Add ligand and parameters to self
    #
    $params->{ligand_file} = $ligand_file; # all ligand inputs are written to file
    $self -> {params} = $params;
    # 
    my @out;
    #
    # Add flexibility to add proteins from PDB selector and via PDB file upload
    #
    if ($params->{protein_input_type} eq 'input_pdb')
        {
        #
        # Stage the PDB data
        #
        my @pdb_list = $self->stage_pdb($params->{input_pdb});

        #
        # And compute
        #

        for my $pdb (@pdb_list)
            {
            my $work = "$work_dir/$pdb->{pdb_id}";
            make_path($work);
            chdir($work);
            $self->compute_pdb($pdb, $ligand_file, $work);
            $self->write_report(\@pdb_list);
            $self->save_output_files($app, $output_dir);
            }
        }
        elsif ($params->{protein_input_type} eq 'user_pdb_file')
        {
            my @ws_pdb_file_list = @{$params->{user_pdb_file}};
            # wrap into a function # stage_ws_pdb
            foreach my $pdb_ws_file (@ws_pdb_file_list) {
                my @parts = split('/', $pdb_ws_file);
                # Get the last element in the array (the file name)
                my $filename = $parts[-1];  
                my @cmd = ("p3-cp", "ws:" . $pdb_ws_file, $self->staging_dir . "/" . basename($filename));
                    print STDERR "Copying the workspace pdb file to stagging dir... @cmd\n";
                    my $ok = IPC::Run::run(\@cmd);
                    if (!$ok)
                        {
                        warn "Error $? copying output with @cmd\n";
                        }
                my $pdb_file = $self->staging_dir . "/" . basename($filename);
                my $protein_id;

                # Open the PDB file for reading
                open (my $fh, '<', $pdb_file) or die "Could not open file '$pdb_file' $!";
                my $protein_id;

                # Try to find the protein ID in the HEADER line
                while (my $line = <$fh>) {
                    # Look for the HEADER line
                    if ($line =~ /^HEADER/) {
                        # Extract the protein ID from the last 4 characters of the line
                        $protein_id = substr($line, 62, 4);
                        $protein_id =~ s/\s+//g;  # Remove any extra whitespace
                        last;  # Exit the loop after finding the HEADER line
                    }
                } 
                close($fh);
                
                # Always use the filename as a fallback if protein ID is not valid
                if (!defined $protein_id || $protein_id eq '') {
                    # Use the filename (without extension) as the fallback protein ID
                    $protein_id = basename($filename, ".pdb");
                    warn "Protein ID not found in HEADER line; using filename '$protein_id' as fallback.";
                }
                # Proceed with the protein ID
                print "Protein ID: $protein_id\n";

                my $work = "$work_dir/$protein_id";
                make_path($work) or die "Could not create directory'$work': $!";

                chdir($work);
                # make pdb object to match the rest of the code
                my $pdb = {
                    "local_path" => $pdb_file,
                    "pdb_id" => $protein_id,
                };
                $self->compute_pdb($pdb, $ligand_file, $work);
                $self->write_report($pdb);
                $self->save_output_files($app, $output_dir);
                }
        }
    return @out;
}


sub write_zero_valid_ligands_report
{
    # my($self, $pdbs) = @_;
    # my($self, $invalid_ligands_file) = @_;
    my($self) = @_;

    my $url_base = $ENV{P3_BASE_URL} // "https://www.bv-brc.org";
    my %vars = (
        # proteins => $pdbs,
		work_dir => $self->{work_dir},
        staging_dir => $self->{staging_dir},
        output_dir => $self->{output_dir},
		ligands => $self->{ligand_name},
		ligand_info => $self->{ligand_info},
        failed_validation => $self->staging_dir . "/invalid_smile_strings.txt",
		params => $self->params,
		output_folder => $self->params->{output_path} . "/." . $self->params->{output_file},
		url_base => $url_base,
		feature_base => "$url_base/view/Feature",
		structure_base => "$url_base/view/ProteinStructure#path",
        bvbrc_logo => "/vol/bvbrc/production/application-backend/bvbrc_docking/bv-brc-header-logo-bg.png"
			);
	# Convert the hash to a JSON string 
	my $json_text = to_json(\%vars, { pretty => 1 });
	#Define the path to the report_data.json file
	my $report_data_path = File::Spec->catfile($self->{work_dir}, "report_data.json");
    my $raw_report_tsv = File::Spec->catfile($self->{work_dir}, "raw_report_data.tsv");
    my $output_html_table = File::Spec->catfile($self->{output_dir}, "docking_results_explorer.html");
    my $output_results_tsv = File::Spec->catfile($self->{output_dir}, "docking_results.tsv");

	# Write the JSON string to the file
	open(my $fh, '>', $report_data_path) or die "Could not open file '$report_data_path': $!";
	print $fh $json_text;
	close($fh);
	print "Analysis data written to $report_data_path\n";

    # Docking report
	my @cmd = (
        "write_docking_html_report",
		"$report_data_path"
	);

    print STDERR "Run: @cmd\n";
    my $ok = IPC::Run::run(\@cmd);
    if (!$ok)
    {
     die "Report command failed $?: @cmd";
    }

    # TSV to HTML viewer
	my @new_cmd = (
		"tsv_to_html",
		"$raw_report_tsv",
        "$output_html_table",
        "$output_results_tsv"
	);
    print STDERR "Run new command: @new_cmd\n";
    my $new_ok = IPC::Run::run(\@new_cmd);
    if (!$new_ok)
    {
     die "HTML Table command failed $?: @new_cmd";
    }
}

sub write_report
{
    my($self, $pdbs) = @_;
    
    my $url_base = $ENV{P3_BASE_URL} // "https://www.bv-brc.org";
    my %vars = (
        proteins => $pdbs,
		work_dir => $self->{work_dir},
        staging_dir => $self->{staging_dir},
        output_dir => $self->{output_dir},
		ligands => $self->{ligand_name},
		ligand_info => $self->{ligand_info},
        failed_validation => $self->{failed_validation},
		results => $self->{result_data},
		params => $self->params,
		output_folder => $self->params->{output_path} . "/." . $self->params->{output_file},
		url_base => $url_base,
		feature_base => "$url_base/view/Feature",
		structure_base => "$url_base/view/ProteinStructure#path",
        bvbrc_logo => "/vol/bvbrc/production/application-backend/bvbrc_docking/bv-brc-header-logo-bg.png"
			);

	# Convert the hash to a JSON string 
	my $json_text = to_json(\%vars, { pretty => 1 });
	#Define the path to the report_data.json file
	my $report_data_path = File::Spec->catfile($self->{work_dir}, "report_data.json");
    my $raw_report_tsv = File::Spec->catfile($self->{work_dir}, "raw_report_data.tsv");
    my $output_html_table = File::Spec->catfile($self->{output_dir}, "docking_results_explorer.html");
    my $output_results_tsv = File::Spec->catfile($self->{output_dir}, "docking_results.tsv");
	
    # Write the JSON string to the file
	open(my $fh, '>', $report_data_path) or die "Could not open file '$report_data_path': $!";
	print $fh $json_text;
	close($fh);
	print "Analysis data written to $report_data_path\n";

	my @cmd = (
		"write_docking_html_report",
		"$report_data_path"
	);

    print STDERR "Run: @cmd\n";
    my $ok = IPC::Run::run(\@cmd);
    if (!$ok)
    {
     die "Report command failed $?: @cmd";
    }

    # TSV to HTML viewer
	my @new_cmd = (
		"tsv_to_html",
		"$raw_report_tsv",
        "$output_html_table",
        "$output_results_tsv"
	);
    print STDERR "Run new command: @new_cmd\n";
    my $new_ok = IPC::Run::run(\@new_cmd);
    if (!$new_ok)
    {
     die "HTML Table command failed $?: @new_cmd";
    }
}

sub compute_pdb
{
    my($self, $pdb, $ligand_file, $work_dir) = @_;
    my $work_out = "$work_dir/out";

    #
    # Determine if our protein is large enough to require a smaller batch_size
    #
    my $residues;
    my $ok = IPC::Run::run(["count-pdb-residues", $pdb->{local_path}], ">", \$residues);
    chomp $residues;
    if (!$ok)
    {
	die "Could not determine residue count for $pdb->{local_path}";
    }
    print STDERR "PDB has $residues residues\n";
    if ($residues > 1024)
    {
	$self->params->{batch_size} = 5 if $self->params->{batch_size} > 5;
	print STDERR "Setting batch_size to " . $self->params->{batch_size} . "\n";
    }
    
    my @batch_size;
    if ($self->params->{batch_size} =~ /\d/)
    {
	@batch_size = ('--batch-size', $self->params->{batch_size});
    }

    #
    # Use simplified diffdock_run.py script
    #
    my @cmd = ('diffdock_run',
	       '--pdb', $pdb->{local_path},
	       '--ligands', $ligand_file,
	       '--outdir', $work_out,
	       @batch_size);

    # Add optional parameters if specified
    if ($self->params->{samples_per_complex}) {
	push @cmd, '--samples-per-complex', $self->params->{samples_per_complex};
    }
    if ($self->params->{inference_steps}) {
	push @cmd, '--inference-steps', $self->params->{inference_steps};
    }
    if ($self->params->{top_n}) {
	push @cmd, '--top-n', $self->params->{top_n};
    }

    print "RUN @cmd\n";
    #
    # Quasi debugging; don't run if the output is already there.
    if (! -d $work_out)
    {
	my $ok = IPC::Run::run(\@cmd);
	
	print "ok=$ok\n";
	if (!$ok)
	{
	    die "Docking run failed\n";
	}
    }

    #
    # Copy output from work dir to the output directory.
    # We preserve ligand output directory from the docking run
    # and just move it into the output directory.
    #
    # While we're there, we load the result.csv and store the contents
    # for later reporting.
    #

    my $out = $self->output_dir . "/$pdb->{pdb_id}";
    make_path($out);

    for my $ligand (@{$self->ligand_name})
    {
	print "Copy $work_out/$ligand to $out\n";
	if (!dircopy("$work_out/$ligand", "$out/$ligand"))
	{
	    die "Failed to copy $work_out/$ligand to $out";
	}
	#
	# Also copy the diffdock logfile.
	#
	copy("$work_out/diffdock_log", "$out/diffdock_log.txt");

	if (-f "$work_out/bad-ligands.txt") {
	    copy("$work_out/bad-ligands.txt", "$out/bad-ligands.txt");
	}
    # File exists even if it is empty - Checking size
    my $staging_dir = $self->staging_dir;
    if (-f "$staging_dir/invalid_smile_strings.txt" && -s "$staging_dir/invalid_smile_strings.txt") {
    copy("$staging_dir/invalid_smile_strings.txt" , "$out/invalid_smile_strings.txt");
    }

	
	my $result_data = csv(in => "$work_out/$ligand/result.csv", headers => 'auto', sep_char => "\t");
	$_->{output_folder} = "$pdb->{pdb_id}/$ligand" foreach @$result_data;
	$self->{result_data}->{$pdb->{pdb_id}}->{$ligand} = $result_data;
    }
}

sub load_ligand_smiles
{
    my($self, $smiles_list) = @_;

    #
    # We just store this in the smiles_list member.
    #
    $self->{smiles_list} = $smiles_list;
    my $staging_dir = $self->staging_dir;
    
    my $file = $self->staging_dir . "/raw_ligands.smi";
    my $validated_ligands_file = $self->staging_dir . "/ligands.smi";

    open(F, ">", $file) or die "Cannot write $file: $!";
    my $row = 0;
    my @new;
    for my $elt_in (@$smiles_list)
    {
	my $elt;
	my $id;
	if (ref($elt_in) eq 'ARRAY')
	{
	    next if @$elt_in == 0;
	    if (@$elt_in == 1)
	    {
		$elt = $elt_in->[0];
	    }
	    else
	    {
		$id = $elt_in->[0];
		$elt = $elt_in->[1];
	    }
	}
	else
	{
	    $elt = $elt_in;
	}
	$id //= sprintf("ligand-%04d", $row + 1);
	print F join("\t", $id, $elt), "\n";
	push(@new, $elt);
	$row++;
    }
    close(F);

    #RUN Validate the input ligands
    my @cmd = (
        "check_input_smile_strings",
        "$staging_dir",
        "$file",
        );
    print STDERR "Run: @cmd\n";
    my $ok = IPC::Run::run(\@cmd);
    if (!$ok)
    {
     die "Input clean up command failed $?: @cmd";
    }
    # Assign the ligand values from the validated ligands.smi
    open(VF, "<", $validated_ligands_file) or die "Cannot read $validated_ligands_file: $!";
    $row = 0;
    @new = ();

    while (my $line = <VF>) {
        chomp $line;
        my ($id, $elt) = split("\t", $line);

        $self->{ligand_map}->{$id} = $row;
        $self->{ligand_name}[$row] = $id;
        $self->{ligand_info}[$row] = { id => $id, idx => $row, smiles => $elt };
        push(@new, $elt);

        $row++;
    }

    $self->{smiles_list} = \@new;
    close(VF);

    if (-e $validated_ligands_file) { # Check if file exists
        if (-s $validated_ligands_file == 0) {  # Check if the file size is zero
            if (-f "$staging_dir/invalid_smile_strings.txt" && -s "$staging_dir/invalid_smile_strings.txt") {
                    # run the report script
                    $self->write_zero_valid_ligands_report();
                    my $output = $self->output_dir;
                    my $workspace_output_path = $self->params->{output_path} . "/." . $self->params->{output_file};
                    my @cmd = ("p3-cp", "--overwrite", "$output/small_molecule_docking_report.html", "ws:" . $workspace_output_path);
                    print STDERR "saving report to workspace... @cmd\n";
                    my $ok = IPC::Run::run(\@cmd);
                    if (!$ok)
                        {
                        warn "Error $? copying output with @cmd\n";
                        }
                    print STDERR "Report uploaded. Exiting before DiffDock runs ";
                    exit 0;
                    }
        } else {
            # Continue and pass the valid ligands
            return $validated_ligands_file;
        }
    } else {
        die "Validated ligands file does not exist.\n";
    }
}


sub load_ligand_library {
    my ($self, $lib_name) = @_;

    #
    # Given a ligand library write ID, Name, and SMILE string to info.txt with names
    # And pass only ID and SMILE string to match the other inputs
    #

    # Open the ligand library or local copy of 3 col ws file in the staging directory
    open(my $fh, '<', $lib_name) or die "Could not open file '$lib_name' $!";
    my $file_contents = do { local $/; <$fh> };
    close($fh);

    my $dat = ($file_contents);
    my $staging_dir = $self->staging_dir;

    # Open info.txt for writing in the staging directory - exsits with ws files and ligand libs
    open(my $info_fh, '>', "$staging_dir/info.txt") or die "Could not open file '$staging_dir/info.txt' $!";

    open(IN, "<", \$dat) or die "Cannot string-open results: $!";
    my @dat;

    while (<IN>) {
        chomp;
        s/^\s*//; # Remove leading whitespace
        next if $_ eq '';
        my @cols = split(/\s+/);

        # Check if there are at least three columns (to accomodate names with spaces without tabs)
        if (scalar @cols >= 3) {
            # Push only ID and SMILES string to @dat
            push(@dat, [$cols[0], $cols[-1]]); # Only pass the ID (first column) and SMILES (last column) to @dat
            # Write all columns to info.txt with tab between only the first and last columns
            print $info_fh $cols[0] . "\t" . join(" ", @cols[1..$#cols-1]) . "\t" . $cols[$#cols] . "\n";
        }
    }
    close($info_fh);

    # Pass the filtered data to load_ligand_smiles
    my $res = $self->load_ligand_smiles(\@dat);
    return $res;
}

sub load_ligand_ws_file
{
    my($self, $ws_file) = @_;

    # Download the file contents as a string from the workspace
    my $dat = $self->app->workspace->download_file_to_string($ws_file);

    # Open the string data for reading
    open(IN, "<", \$dat) or die "Cannot string-open results: $!";
    my @dat;
    my $columns = 0;

    # Read the file and determine the number of columns
    while (<IN>)
    {
        chomp;
        s/^\s*//;
        next if $_ eq '';

        # Split by tabs or multiple spaces
        my @cols = split(/\t|\s{2,}/);
        $columns = scalar(@cols) if $columns == 0;
        push(@dat, [@cols]);
    }
    close(IN);

    if ($columns == 2) {
        # If there are two columns, pass the data to load_ligand_smiles
        return $self->load_ligand_smiles(\@dat);
    }
    elsif ($columns >= 3) {
        # If there are three columns, write first and third columns to a local file

        # Define the local file path in the staging directory
        my $local_file = $self->staging_dir . "/three_col_ws_file.txt";

        # Open the local file for writing
        open(my $local_fh, '>', $local_file) or die "Could not open file '$local_file' $!";

        # Write data to the local file to give ligand library
        foreach my $row (@dat) {
            print $local_fh join("\t", @$row) . "\n";  # Use join to print all columns in the row
        }
        close($local_fh);

        # Pass the local file path to load_ligand_library like a library
        return $self->load_ligand_library($local_file);
    }
    else {
        die "Unexpected number of columns ($columns) in file.";
    }
}


sub stage_pdb
{
    my($self, $pdb_list) = @_;
    
    my $qry = join(",", @$pdb_list);
    my @res = $self->api->query('protein_structure', ['in', 'pdb_id', '(' . $qry . ')'], ['select', 'pdb_id,gene,product,method,patric_id,title,file_path']);

    print Dumper(@res);
    my %res = map { ($_->{pdb_id} => $_ ) } @res;

    my @out;
    for my $pdb (@$pdb_list)
    {
	my $ent = $res{$pdb};
	if (!$ent)
	{
	    die "No pdb found for $pdb\n";
	}
	my $path = "/vol/structure/$ent->{file_path}";
	if (! -s $path)
	{
	    die "Cannot find $path\n";
	}
	my $local = $self->staging_dir . "/" . basename($path);
	copy($path, $local);
	$ent->{local_path} = $local;
	$self->{pdb_info}->{$pdb} = $ent;
	push(@out, $ent);
    }
    return @out;
}


sub preflight
{
    my ($app, $app_def, $raw_params, $params) = @_;

    my $resource;
    
    my %resource_map = (
        named_library => {
            "experimental_drugs" => { mem => "128G", time => 60 * 60 * 18 }, # 18 hours
            "approved-drugs" => { mem => "128G", time => 60 * 60 * 16 }, # 16 hours
            "small_db" => { mem => "128G", time => 60 * 60 * 1  }, # 1 hours
        },
        smiles_list => { mem => "128G", time => 60 * 60 * 4 }, # 4 hours
        ws_file     => { mem => "128G", time => 60 * 60 * 4 }, # 4 hours
    );

    if (exists $params->{ligand_library_type}) {
        if ($params->{ligand_library_type} eq 'named_library' && exists $params->{ligand_named_library}) {
            $resource = $resource_map{named_library}{$params->{ligand_named_library}};
        } elsif ($params->{ligand_library_type} eq 'smiles_list') {
            $resource = $resource_map{smiles_list};
        } elsif ($params->{ligand_library_type} eq 'ws_file') {
            $resource = $resource_map{ws_file};
        } else {
            die "Unknown ligand library type selected: $params->{ligand_library_type}";
        }
    } else {
        die "Ligand library type not specified";
    }

    my $mem = $resource->{mem};
    my $runtime = $resource->{time};

    my $pf = {
        cpu => 8,
        memory => $mem,
        runtime => $runtime,
        policy_data => { gpu_count => 1, partition => 'gpu', constraint => 'V100' },
    };
    return $pf;
}

sub save_output_files
{
    my($self, $app, $output) = @_;
    my %suffix_map = (
		      bai => 'bai',
		      bam => 'bam',
		      #
		      # The result.csv generated by the app is really tab separated.
		      #
		      csv => 'tsv',
		      depths => 'txt',
		      err => 'txt',
		      fasta => "contigs",
		      html => 'html',
		      out => 'txt',
		      txt => 'txt',
		      png => 'png',
		      pdb => 'pdb',
		      tsv => 'tsv',);

    my @suffix_map = map { ("--map-suffix", "$_=$suffix_map{$_}") } keys %suffix_map;

    if (opendir(D, $output))
    {
	while (my $p = readdir(D))
	{
	    next if ($p =~ /^\./);
	    my @cmd = ("p3-cp", "--overwrite", "--recursive", @suffix_map, "$output/$p", "ws:" . $app->result_folder);
	    print STDERR "saving files to workspace... @cmd\n";
	    my $ok = IPC::Run::run(\@cmd);
	    if (!$ok)
	    {
		warn "Error $? copying output with @cmd\n";
	    }
	}
    closedir(D);
    }
}

1;
