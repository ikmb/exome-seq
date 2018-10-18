file = ARGV.shift
lines = `zcat #{file}`.split("\n")

lines.each do |line|
	e = line.strip.split(/\s+/)
	puts "#{e[0]}:#{e[1]}-#{e[2]}"
end
