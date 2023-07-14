#!/usr/bin/env ruby
# == NAME
# vep2xlsx.rb
#
# == AUTHOR
#  Marc Hoeppner, mphoeppner@gmail.com

require 'optparse'
require 'ostruct'
require 'fast_excel'

### Define modules and classes here

class VEPEntry

    attr_accessor :seq, :pos, :rsid, :ref, :alt, :qual, :filter, :info, :format_info, :samples, :csq

    def initialize(line,csq_header)

        elements = line.strip.split("\t")

        @seq,@pos,@rsid,@ref,@alt,@qual,filter,@info,@format_info = elements[0..8]
        @samples = elements[9..-1]
        @filter = filter.split(";")[0]

        @seq.gsub!(/^chr/, '')
        @csq = []

        csq_field = @info.split(";").find { |e| e.include?("CSQ") }
        if csq_field
            tmp = csq_field.split("=")[-1]

            csq = tmp.split(",").each do |c_e|

                data = c_e.split("|")
                this_csq = {}

                csq_header.each_with_index do |c,i|

                    this_data = data[i]
                    if this_data && this_data.match(/^[0-9]*\.[0-9]*/)
                        this_data.gsub!(".",",")
                    end
                    this_csq[c] = this_data

                end

                @csq << this_csq
        
            end

        end

    end

	def genotypes

		answer = []
		alleles = [ @ref ]
        @alt.split(",").each { |a| alleles << a}

		genotype_column_index = @format_info.split(":").index("GT")

		self.samples.each do |s|
			genotype_data = s.split(":")[genotype_column_index].split("/").collect {|gt| gt.to_i }
			
            calls = []
            genotype_data.each do |gt|
                calls << alleles[gt]
            end
            answer << calls.join("/")
		end

        return answer
	end

    def genotype(idx)

        genotype_column_index = @format_info.split(":").index("GT")

        return [ ] if @samples[idx].split(":")[genotype_column_index] == "./."

        genotype_data = @samples[idx].split(":")[genotype_column_index].split("/").collect{|gt| gt.to_i}

        alleles = [ @ref ]
        @alt.split(",").each { |a| alleles << a}

        calls = []
        genotype_data.each do |gt|
            calls << alleles[gt]
        end

        return calls

    end

    def annotations(idx)
        results = {}
        keys = self.format_info.split(":")
        self.samples[idx].split(":").each_with_index do |f,i|
            results[keys[i]] = f
        end
        return results
    end

    def genotype_depth

        depth_column_index = @format_info.split(":").index("DP")
        return @sample.split(":")[depth_column_index].to_i
    end

end


### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.on("-i","--infile", "=INFILE","Infile file") {|argument| options.infile = argument }
opts.on("-o","--outfile", "=OUTFILE","Output file") {|argument| options.outfile = argument }
opts.on("-h","--help","Display the usage information") {
    puts opts
    exit
}

opts.parse! 

options.infile ? input_stream = File.open(options.infile,"r") : input_stream = $stdin

color = {
	"even" => "EBEBEB",
	"uneven" => "D6D6D6"
}

workbook = FastExcel.open(options.outfile, constant_memory: true)

workbook.default_format.set(
    font_size: 10, # user's default
    font_family: "Arial"
)

bold                = workbook.bold_format
sheet               = workbook.add_worksheet("VEP results")
sheet.auto_width    = true

even                = workbook.add_format(bg_color: color["even"])
uneven              = workbook.add_format(bg_color: color["uneven"])

csq_header          = []
header              = []

counter             = 0

while (line = input_stream.gets)

    line.strip!

    # The line defining the VEP data structure
    if  line.match(/^##INFO\=\<ID\=CSQ.*/) 
        line.split(" ")[-1].split("|").each do |e|
            csq_header << e.strip.split(/"/)[0]
        end            
    end

    # The VCF column header
    if line.match(/^#CHR.*/)
        header = line.gsub(/^#/, '').split("\t")[0..8]
        header_samples = line.gsub(/^#/, '').split("\t")[9..-1]

        full_header = header + csq_header + header_samples
        full_header = full_header.select {|h| h != "FORMAT" && h != "INFO" }

        sheet.append_row(full_header,bold)
		
    end    
    
    # Line is a comment, skip
    next if line.match(/^#.*/)

    e = VEPEntry.new(line.strip,csq_header)

    # :seq, :pos, :rsid, :ref, :alt, :qual, :filter, :info, :format_info, :sample, :csq
    samples = e.samples

    # add the basic variant information

    e.csq.each do |csq_t|

        values = [ e.seq , e.pos, e.rsid.gsub(";",","), e.ref, e.alt, e.qual, e.filter ]

        counter += 1

        csq_header.each do |cq|
            values << csq_t[cq]
        end

        values = values + e.genotypes

        counter.even? ? bg = even : bg = uneven

        sheet.append_row(values,bg)

    end
    
end

workbook.close
