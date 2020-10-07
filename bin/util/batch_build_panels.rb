this_dir = __dir__

container = "shub://ikmb/ensembl-api:99"

assemblies = [ "hg19", "GRCh37", "GRCh38"  ]

root_folder = "#{this_dir}/../../assets/panels"


lookup = { "hg19" => "hg19", "GRCh37" => "b37", "GRCh38" => "hg38" }

panel_lists = Dir["#{this_dir}/../../assets/panels/gene_lists/*.txt"]

perl_script = Dir["#{this_dir}/../ensembl_panel2bed.pl"]

system("mkdir -p $HOME/exome_assets_backup")

assemblies.each do |assembly|

	backup_command = "cp -R #{root_folder}/#{assembly} $HOME/exome_assets_backup/"
	warn backup_command	
	warn "Generating lists for assembly version #{assembly}"

	panel_lists.each do |pl|
		name = pl.split("/")[-1].gsub(".txt", ".#{assembly}.bed")
		command = "singularity exec -B /work_ifs #{container} perl #{perl_script} --assembly #{assembly} --list #{pl} > #{name}"

		warn command
	end
end
