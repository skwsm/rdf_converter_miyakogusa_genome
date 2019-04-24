#!/usr/bin/env ruby

print "@prefix m2r: <http://med2rdf.org/ontology/med2rdf#> .\n"
print "@prefix faldo: <http://biohackathon.org/resource/faldo#> .\n"
print "@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .\n"
print "@prefix dcterms: <http://purl.org/dc/terms/> .\n"
print "@prefix obo: <http://purl.obolibrary.org/obo/> .\n"
print "@prefix so: <http://purl.obolibrary.org/obo/so#> .\n"
print "@prefix sio: <http://semanticscience.org/resource/> .\n"
print "@prefix miya: <http://www.kazusa.or.jp/lotus/predict/cgi-bin/orfinfo.cgi?db=Lj3.0&id=> .\n"
print "\n"

f_gff = open(ARGV.shift) # GFF3


class GFF2RDF
  def initialize()
    @gene2transcripts = {}
    @transcript2exons = {}
    @transcript2cdss   = {}
  end
  attr_accessor :gene2transcripts, :transcript2exons, :transcript2cdss

  def parse_annotation(annot)
    annot_hash = {}
    ary = annot.split(";")
    ary.each do |e|
      case e
      when /ID\=(\S+)$/
        annot_hash[:identifier] = $1
      when /Parent\=(\S+)$/
        annot_hash[:parent] = $1
      when /Name=(\S+)$/
        annot_hash[:name] = $1
      when /Type=(\S+)$/
        annot_hash[:type] = $1
      else
      end
    end
    annot_hash
  end

  def construct_gene_models(gff_line)
    annot_hash = parse_annotation(gff_line[:annotation])
    gff_line[:identifier] = annot_hash[:identifier] if annot_hash.key?(:identifier)
    gff_line[:name]       = annot_hash[:name]       if annot_hash.key?(:name)
    gff_line[:parent]     = annot_hash[:parent]     if annot_hash.key?(:parent)

    case gff_line[:feature]
    when "gene"
      unless @gene2transcripts.key?(annot_hash[:identifier])
        @gene2transcripts[annot_hash[:identifier]] = gff_line
      end
    when "mRNA"
      if @gene2transcripts.key?(annot_hash[:parent])
        if @gene2transcripts[annot_hash[:parent]].key?(:transcripts)
          @gene2transcripts[annot_hash[:parent]][:transcripts] << gff_line
        else
          @gene2transcripts[annot_hash[:parent]][:transcripts] = [gff_line]
        end
      end
    when "CDS"
      if @transcript2cdss.key?(annot_hash[:parent])
        @transcript2cdss[annot_hash[:parent]] << gff_line
      else
        @transcript2cdss[annot_hash[:parent]] = [gff_line]
      end
    when "exon"
      if @transcript2exons.key?(annot_hash[:parent])
        @transcript2exons[annot_hash[:parent]] << gff_line
      else
        @transcript2exons[annot_hash[:parent]] = [gff_line]
      end
    else
    end
  end

  def print_rdf
    @gene2transcripts.each do |gene, gene_hash|
      print "miya:#{gene} a m2r:Gene ;\n"
      case gene_hash[:feature]
      when "protein_coding"
        print "  a obo:SO_0001217 ;\n"
      when "processed_transcript"
        print "  a obo:SO_0001503 ;\n"
      when "rRNA"
        print "  a obo:SO_0001637 ;\n"
      when "tRNA"
        print "  a obo:SO_0001272 ;\n"
      end
      print "  dcterms:identifier \"#{gene}\" ;\n"
      print "  rdfs:label \"#{gene_hash[:name]}\" ;\n"
      print "  faldo:location [\n"
      print "    a faldo:Region ;\n"
      print "    faldo:begin [\n"
      print "      a faldo:ExactPosition ;\n"
      print "      a faldo:ForwardStrandPosition ;\n" if gene_hash[:direct] == 1
      print "      a faldo:ReverseStrandPosition ;\n" if gene_hash[:direct] == -1
      print "      faldo:position #{gene_hash[:start]} ;\n"
      print "      faldo:reference miya:#{gene_hash[:chr_num]}\n"
      print "    ] ;\n"
      print "    faldo:end [\n"
      print "      a faldo:ExactPosition ;\n"
      print "      a faldo:ForwardStrandPosition ;\n" if gene_hash[:direct] == 1
      print "      a faldo:ReverseStrandPosition ;\n" if gene_hash[:direct] == -1
      print "      faldo:position #{gene_hash[:stop]} ;\n"
      print "      faldo:reference miya:#{gene_hash[:chr_num]}\n"
      print "    ]\n"
      print "  ] .\n"
      gene_hash[:transcripts].each do |transcript|
        print "miya:#{transcript[:identifier]}\n"
        print "  dcterms:identifier \"#{transcript[:identifier]}\" ;\n"
        print "  so:transcribed_from miya:#{transcript[:parent]} ;\n"
        print "  faldo:location [\n"
        print "    a faldo:Region ;\n"
        print "    faldo:begin [\n"
        print "      a faldo:ExactPosition ;\n"
        print "      a faldo:ForwardStrandPosition ;\n" if transcript[:direct] == "+"
        print "      a faldo:ReverseStrandPosition ;\n" if transcript[:direct] == "-"
        print "      faldo:position #{transcript[:start]} ;\n"
        print "      faldo:reference miya:#{transcript[:chr_num]}\n"
        print "    ] ;\n"
        print "    faldo:end [\n"
        print "      a faldo:ExactPosition ;\n"
        print "      a faldo:ForwardStrandPosition ;\n" if transcript[:direct] == "+"
        print "      a faldo:ReverseStrandPosition ;\n" if transcript[:direct] == "-"
        print "      faldo:position #{transcript[:start]} ;\n"
        print "      faldo:reference miya:#{transcript[:end]}\n"
        print "    ]\n"
        print "  ] .\n"
        print "\n"
        @transcript2exons[transcript[:identifier]].each do |exon|
          exon_id = exon[:identifier]
           /[n\.\:](\d+)$/ =~ exon_id
          exon_num = $1
          print "miya:#{transcript[:identifier]} so:has_part miya:#{exon_id} ;\n"
          print "  sio:SIO_000974 miya:#{transcript[:identifier]}\\#Exon_#{exon_num} .\n\n"
          print "miya:#{exon_id} a obo:SO_0000147 ;\n"
          print "  dcterms:identifier \"#{exon_id}\" ;\n"
          print "  so:part_of miya:#{transcript[:identifier]} ;\n"
          print "  faldo:location [\n"
          print "    a faldo:Region ;\n"
          print "    faldo:begin [\n"
          print "      a faldo:ExactPosition ;\n"
          print "      a faldo:ForwardStrandPosition ;\n" if exon[:direct] == "+"
          print "      a faldo:ReverseStrandPosition ;\n" if exon[:direct] == "-"
          print "      faldo:position #{exon[:start]} ;\n"
          print "      faldo:reference miya:#{exon[:chr_num]}\n"
          print "    ] ;\n"
          print "    faldo:end [\n"
          print "      a faldo:ExactPosition ;\n"
          print "      a faldo:ForwardStrandPosition ;\n" if exon[:direct] == "+"
          print "      a faldo:ReverseStrandPosition ;\n" if exon[:direct] == "-"
          print "      faldo:position #{exon[:start]} ;\n"
          print "      faldo:reference miya:#{exon[:chr_num]}\n"
          print "    ]\n"
          print "  ] .\n"
          print "miya:#{transcript[:identifier]}\\#Exon_#{exon_num} \n"
          print "  sio:SIO_000300 #{exon_num} ;\n"
          print "  sio:SIO_000628 miya:#{transcript[:identifier]} ;\n"
          print "  a sio:SIO_001261 .\n"
          print "\n"
        end
        if @transcript2cdss.key?(transcript[:identifier])
        @transcript2cdss[transcript[:identifier]].each do |cds|
          cds_id = cds[:identifier]
          /[Ss\:\.](\d+)$/ =~ cds_id
          cds_num = $1
          print "miya:#{transcript[:identifier]} so:has_part miya:#{cds_id} ;\n"
          print "  sio:SIO_000974 miya:#{transcript[:identifier]}\\#CDS_#{cds_num} .\n\n"
          print "miya:#{cds_id} a obo:SO_0000316 ;\n"
          print "  dcterms:identifier \"#{cds_id}\" ;\n"
          print "  so:part_of miya:#{transcript[:identifier]} ;\n"
          print "  miya:frame #{cds[:frame]} ;\n"
          print "  faldo:location [\n"
          print "    a faldo:Region ;\n"
          print "    faldo:begin [\n"
          print "      a faldo:ExactPosition ;\n"
          print "      a faldo:ForwardStrandPosition ;\n" if cds[:direct] == "+"
          print "      a faldo:ReverseStrandPosition ;\n" if cds[:direct] == "-"
          print "      faldo:position #{cds[:start]} ;\n"
          print "      faldo:reference miya:#{cds[:chr_num]}\n"
          print "    ] ;\n"
          print "    faldo:end [\n"
          print "      a faldo:ExactPosition ;\n"
          print "      a faldo:ForwardStrandPosition ;\n" if cds[:direct] == "+"
          print "      a faldo:ReverseStrandPosition ;\n" if cds[:direct] == "-"
          print "      faldo:position #{cds[:start]} ;\n"
          print "      faldo:reference miya:#{cds[:chr_num]}\n"
          print "    ]\n"
          print "  ] .\n"
          print "miya:#{transcript[:identifier]}\\#CDS_#{cds_num} \n"
          print "  sio:SIO_000300 #{cds_num} ;\n"
          print "  sio:SIO_000628 miya:#{transcript[:identifier]} ;\n"
          print "  a sio:SIO_001261 .\n"
          print "\n"
        end
        end
      end
    end
  end
end

gm = GFF2RDF.new

while line = f_gff.gets
  keys = [:chr_num, :gene_model_type, :feature, :start, :stop, :score, :direct, :frame, :annotation]
  vals = line.chomp.split("\t")

  gff_line = Hash[keys.zip(vals)]
  gm.construct_gene_models(gff_line)
end

gm.print_rdf 


