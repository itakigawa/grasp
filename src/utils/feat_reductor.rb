#!/usr/bin/ruby
# ruby feat_reductor.rb output_smol05/feat.instances output_smol05/pair.patterns 10000 > nonredundant_10000.patterns

f = open(ARGV[0])
g = open(ARGV[1])
num = ARGV[2].to_i

h = Hash.new
c = 0
removed = 0

while line = f.gets
  rank,signature = line.chomp.split("|")
  if not h.has_key?(signature)
    h[signature] = rank.to_i
    c += 1
  else
    removed += 1
    #STDERR.puts "#{rank} removed because of #{h[signature]}"
  end
  if c == num
    break
  end
end

STDERR.puts "#{removed} removed."

export = h.values.sort

l = 0
a = export.shift
while line = g.gets
  l += 1
  if l == a
    puts line.chomp
    a = export.shift
    if a == nil
      break
    end
  end
end

