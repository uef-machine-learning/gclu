#!/usr/bin/ruby

require "tempfile"
require "ripl"

def svg_head(maxx,maxy)

# <rect x="0" y="0" width="#{(maxx*1.2).to_i}" height="#{(maxy*1.2).to_i}" fill="none"/>
  return <<-SVG00
<?xml version="1.0" encoding="utf-8"  standalone="no"?>
<svg 
 width="#{maxx}" height="#{maxy}"
 viewBox="0 0 #{maxx} #{maxy}"
 xmlns="http://www.w3.org/2000/svg"
 xmlns:xlink="http://www.w3.org/1999/xlink"
>

<g>

<defs>

    <circle id='gpDot' r='0.5' stroke-width='0.5'/>
    <path id='gpPt0' stroke-width='0.267' stroke='currentColor' d='M-1,0 h2 M0,-1 v2'/>
    <path id='gpPt1' stroke-width='0.267' stroke='currentColor' d='M-1,-1 L1,1 M1,-1 L-1,1'/>
    <path id='gpPt2' stroke-width='0.267' stroke='currentColor' d='M-1,0 L1,0 M0,-1 L0,1 M-1,-1 L1,1 M-1,1 L1,-1'/>
    <rect id='gpPt3' stroke-width='0.267' stroke='currentColor' x='-1' y='-1' width='2' height='2'/>
    <rect id='gpPt4' stroke-width='0.267' stroke='currentColor' fill='currentColor' x='-1' y='-1' width='2' height='2'/>
    <circle id='gpPt5' stroke-width='0.267' stroke='currentColor' cx='0' cy='0' r='1'/>
    <circle id='gpPt15' stroke-width='0.267' stroke='currentColor' cx='0' cy='0' r='0.5'/>
    <use xlink:href='#gpPt15' id='gpPt16' fill='currentColor' stroke='none'/>
    <use xlink:href='#gpPt5' id='gpPt6' fill='currentColor' stroke='none'/>
    <path id='gpPt7' stroke-width='0.267' stroke='currentColor' d='M0,-1.33 L-1.33,0.67 L1.33,0.67 z'/>
    <use xlink:href='#gpPt7' id='gpPt8' fill='currentColor' stroke='none'/>
    <use xlink:href='#gpPt7' id='gpPt9' stroke='currentColor' transform='rotate(180)'/>
    <use xlink:href='#gpPt9' id='gpPt10' fill='currentColor' stroke='none'/>
    <use xlink:href='#gpPt3' id='gpPt11' stroke='currentColor' transform='rotate(45)'/>
    <use xlink:href='#gpPt11' id='gpPt12' fill='currentColor' stroke='none'/>
    <path id='gpPt13' stroke-width='0.267' stroke='currentColor' d='M0,1.330 L1.265,0.411 L0.782,-1.067 L-0.782,-1.076 L-1.265,0.411 z'/>
    <use xlink:href='#gpPt13' id='gpPt14' fill='currentColor' stroke='none'/>
    <filter id='textbox' filterUnits='objectBoundingBox' x='0' y='0' height='1' width='1'>
      <feFlood flood-color='white' flood-opacity='1' result='bgnd'/>
      <feComposite in='SourceGraphic' in2='bgnd' operator='atop'/>
    </filter>
    <filter id='greybox' filterUnits='objectBoundingBox' x='0' y='0' height='1' width='1'>
      <feFlood flood-color='lightgrey' flood-opacity='1' result='grey'/>
      <feComposite in='SourceGraphic' in2='grey' operator='atop'/>
    </filter>
</defs>
SVG00
end

def svg_end()
  return <<-SVG00
</svg>
SVG00
end

def svg_point(x, y, type, color)
  opacity = "0.50"
  scale = 7.45
  return "<use xlink:href='#{type}' transform='translate(#{x},#{y}) scale(#{scale})' color='#{color}' opacity='#{opacity}'/>\n"
end

def read_pa(pafn)
  pa = File.readlines(pafn).collect { |x| x.strip }

  if pa[0].match(/VQ PARTITIONING/)
    header_end = -1
    (1..10).each { |i| header_end = i if /\----+/.match(pa[i]) }
    return nil if header_end == -1
    labels = pa[(header_end + 1)..-1].collect { |x| x.to_i }
    k = pa[1].to_i
  else
    labels = pa.collect { |x| x.to_i }
    k = pa.uniq.size
  end
  # Ripl.start :binding => binding

  return [labels, k]
end

def read_vectors(infn, limit = 0)
  arr = []
  infn = File.expand_path(infn)
  File.foreach(infn).with_index do |line, line_num|
    break if limit > 0 and line_num > limit
    # puts "#{line_num}: #{line}"
    a = line.strip.split(/\s+/)
    b = a.collect { |x| x.to_f }
    arr << b
  end
  return arr
end

def create_clu_viz(dsfn, partfn, outfn)
  seed = rand(10000)
  seed = 3462

  (part, k) = read_pa(partfn)
  # puts "k=#{k}"
  ds = read_vectors(dsfn)
  maxx = ds.collect{|x|x[0]}.max
  maxy = ds.collect{|x|x[1]}.max
  
  scalef = 50.0
  scalef = 0.001
  scalef = 100.0/maxy
  scalef = 0.002
  scalef = 1.0
  scalef = 1000.0/[maxx,maxy].max


  # tmpf = Tempfile.new(["clu", ".svg"])
  # Ripl.start :binding => binding

  outf = File.open(outfn, "w")

  # outf.write svg_head(maxx*scalef,maxy*scalef)
  outf.write svg_head(1000,1000)
  point_types = {}
  brush = ["#gpPt6", "#gpPt4", "#gpPt8"]

  # seed=2079
  srand(seed)
  # puts "seed=#{seed}"

  brush = ["#gpPt6", "#gpPt4", "#gpPt8", "#gpPt0", "#gpPt1", "#gpPt2", "#gpPt16"].shuffle

  # brush = ["#gpPt0"] # Cross
  # brush = ["#gpPt1"] #cross diagonal
  # brush = ["#gpPt2"] # Star
  # brush = ["#gpPt16"] # Dot

  colstr = "215.3816  185.9709  209.9076
  125.4288   39.8879  123.6632
   55.8716  198.9620   39.3525
  126.7435  213.6287  205.7354
  200.9577  195.0321   37.5021
   40.5435  101.1956  124.2205
  201.9148   48.6501   47.1747
  183.3573  212.0293  119.5509
  127.7739  116.1867  201.2917
   59.5821  207.5048  126.4313
  210.8995  114.1225  121.1687
   53.4909   45.8396  206.5980
  203.3110   46.7799  203.0992
   53.6669   49.7683   44.5671
   38.8968  182.8289  212.9211
  123.8039  128.0467   54.6903"

  col_lines = colstr.lines.shuffle
  if k > 15
    10.times { |x|
      col_lines << 3.times.collect { |x| r = rand(255) }.join(" ")
    }
  end

  colors = col_lines.collect { |x| a = x.split(" ").collect { |y| (y.to_i * 0.8).to_i }.join(","); "rgb(#{a})" }
  # puts colors
  colors.shuffle!

  for i in 1..(k)
    # pt = brush[i % brush.size]
    pt = brush[rand(brush.size)]

    cmax = 255
    intensity = rand(0.4..1.0)

    colvec = 3.times.collect { |x| a = rand(cmax); r = 0 if rand() < 0.2; r = 255 if rand() < 0.2; r = a * intensity; r.to_i }
    col = "rgb(" + colvec.join(",") + ")"
    col = colors[i - 1]
    # col = "rgb(#{col,#{rand(cmax)},#{rand(cmax)})"
    # col = "rgb(#{rand(cmax)},#{rand(cmax)},#{rand(cmax)})"
    point_types[i] = [pt, col]
  end
  # puts point_types.inspect

  types = ["triangle", "rectangle", "circle", "o", "triangle2", "o2"]
  outf.write '<g fill="none" color="black" stroke="currentColor" stroke-width="1.00" stroke-linecap="butt" stroke-linejoin="miter">'
  dswp = ds.collect.with_index { |x, i| a = x.clone; a << part[i] }
  dswp = dswp.shuffle

  maxy = ds.collect { |x| x[1] }.max
  for (x, y, part_id) in dswp
    pt = point_types[part_id][0]
    col = point_types[part_id][1]
    # outf.write svg_point(scalef * x, scalef * (maxy-y), t)
    outf.write svg_point(scalef * x, scalef * (maxy - y), pt, col)
  end
  outf.write "</g>"
  outf.write "</g>"

  outf.write svg_end
  outf.close
end

dsfn = ARGV[0]
partfn = ARGV[1]
outfn = ARGV[2]
create_clu_viz(dsfn, partfn, outfn)
