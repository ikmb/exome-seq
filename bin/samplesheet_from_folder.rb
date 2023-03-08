#!/bin/env ruby

require 'optparse'
require 'ostruct'

### Define modules and classes here

### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.banner = "Reads Fastq files from a folder and writes a sample sheet to STDOUT"
opts.separator ""
opts.on("-f","--folder", "=FOLDER","Folder to scan") {|argument| options.folder = argument }
opts.on("-c","--centre", "=CENTRE","Name of sequencing centre") {|argument| options.centre = argument }
opts.on("-p","--platform", "=PLATFORM","Name of sequencing instrument") {|argument| options.platform = argument }
opts.on("-s","--sanity", "Perform sanity check of md5 sums") { options.sanity = true }
opts.on("-l","--lookup", "=LOOKUP", "Lookup file with lib id <> other id") {|argument| options.lookup = argument }
opts.on("-h","--help","Display the usage information") {
 puts opts
 exit
}

opts.parse! 

abort "Folder not found (#{options.folder})" unless File.directory?(options.folder)

date = Time.now.strftime("%Y-%m-%d")
options.centre ? center = options.centre : center = "IKMB"

lookup = {}
if options.lookup
	IO.readlines(options.lookup).each do |line|
		key,value = line.strip.split("\t")
		lookup[key] = value
	end
end

fastq_files = Dir["#{options.folder}/*_R*.fastq.gz"]

groups = fastq_files.group_by{|f| f.split("/")[-1].split(/_/)[0..1].join("_") }

warn "Building input sample sheet from FASTQ folder"
warn "Performing sanity check on md5sums" if options.sanity

options.platform ? sequencer = options.platform : sequencer = "NovaSeq6000"

puts "patient;sample;status;library;readgroup;platform_unit;center;date;R1;R2"

individuals = []
samples = []

# group = the library id, may be split across lanes
groups.each do |group, files|

	warn "...processing library #{group}"

	pairs = files.group_by{|f| f.split("/")[-1].split(/_R[1,2]/)[0] }

	pairs.each do |p,reads|

        	left,right = reads.sort.collect{|f| File.absolute_path(f)}
	
		abort "This sample seems to not be a set of PE files! #{p}" unless left && right

		if options.sanity
			Dir.chdir(options.folder) {
				[left,right].each do |fastq|
					fastq_simple = fastq.split("/")[-1].strip
					raise "Aborting - no md5sum found for fastq file #{fastq}" unless File.exists?(fastq_simple + ".md5")
					status = `md5sum -c #{fastq_simple}.md5`
					raise "Aborting - failed md5sum check for #{fastq}" unless status.strip.include?("OK")
				end
			}
		end

		# H26247-L3_S1_L001_R1_001_fastqc.html
		lims_id = group.split("_")[0]
        	library = group.split("_")[1]
        	sample = library
		individual = library

		if lookup.has_key?(library)
			sample = "#{lookup[library]}"
		end

		individuals << individual
		samples << sample

        	e = `zcat #{left} | head -n1 `
		header = e

        	instrument,run_id,flowcell_id,lane,tile,x,y = header.split(" ")[0].split(":")

		index = header.split(" ")[-1].split(":")[-1]
        	readgroup = flowcell_id + "." + lane + "." + library 

        	pgu = flowcell_id + "." + lane + "." + index

        	puts "#{individual};#{sample};0;#{library};#{readgroup};#{pgu};#{center};#{date};#{left};#{right}"
	end
end


warn "Found: #{individuals.uniq.length} Patients."
warn "Found: #{samples.uniq.length} Samples."
warn "If these numbers do not seem right, please re-check the file naming and manually fix the samplesheet."

