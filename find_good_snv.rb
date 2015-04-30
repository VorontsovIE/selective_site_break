require_relative 'perfectosape/results'
require_relative 'sequence_with_snv'
require_relative 'site_list'
require_relative 'output_configuration'

require 'tempfile'
require 'optparse'
require 'set'

def with_temp_file(filename, &block)
  temp_file = Tempfile.new(filename)
  yield temp_file
ensure
  temp_file.close
  temp_file.unlink
end

def make_seq_with_snv_by_position(sequence, pos, mut_nucl)
  raise 'Bad mutation nucleotide'  if mut_nucl.to_s.upcase == sequence[pos].to_s.upcase
  SequenceWithSNV.new(sequence[0,pos], [sequence[pos], mut_nucl], sequence[(pos + 1)..-1])
end

def sequence_with_snv_by_substitution_name(sequence, substitution_name)
  pos, mut_nucleotide = substitution_name.chomp.split(',')
  make_seq_with_snv_by_position(sequence, pos.to_i, mut_nucleotide)
end

def possible_mutation_positions(sequence, ignore_positions: [])
  sequence.length.times.reject{|pos|
    ignore_positions.include?(pos)
  }.flat_map{|pos|
    orig_nucl = sequence[pos]
    mut_nucleotides = ['A','C','G','T'] - [orig_nucl]
    mut_nucleotides.map{|mut_nucl|
      substitution_name = "#{pos},#{mut_nucl}"
      seq_w_snv = make_seq_with_snv_by_position(sequence, pos, mut_nucl)
      [substitution_name, seq_w_snv]
    }
  }.to_h
end

def mutated_sites_by_substitution_name(sequence, pvalue_cutoff: 0.0005, ignore_positions: [])
  with_temp_file 'seq_generated.txt' do |f|
    possible_mutation_positions(sequence, ignore_positions: ignore_positions).each{|substitution_name, seq_w_snv|
      f.puts "#{substitution_name}\t#{seq_w_snv}"
    }
    f.close

    cmd = "java -cp ape.jar ru.autosome.perfectosape.SNPScan ./motif_collection/ #{f.path} --precalc ./motif_collection_thresholds --fold-change-cutoff 1 --pvalue-cutoff 1"
    IO.popen(cmd) do |pipe|
      PerfectosAPE::Result.each_in_stream(pipe).to_a.group_by(&:variant_id)
    end
  end
end

def select_sites_ovelapping_snv(sites_by_substitution, snv_pos)
  sites_by_substitution.map{|variant_id, sites|
    substitution_pos = variant_id.split(",").first.to_i
    selected_sites = sites.select{|site|
      actual_site_position_start = substitution_pos + site.pos_1
      actual_site_position_end = actual_site_position_start + site.length
      (actual_site_position_start...actual_site_position_end).include?(snv_pos)
    }

    [variant_id, selected_sites]
  }.to_h
end

# pattern is space-separated
def choose_motif_names_by_pattern(motif_names, motif_patterns)
  patterns = motif_patterns.map{|pattern| Regexp.new(pattern) }
  motif_names.select{|motif_name|
    patterns.any?{|pattern| pattern.match(motif_name) }
  }.to_set
end


def output_summary_of_sites(sequence, sites_by_substitution,
                            fold_change_cutoff: , pvalue_cutoff: ,
                            to_be_disrupted: , to_be_preserved: ,
                            output_configuration:, stream: $stdout
                          )
  sites_by_substitution.map{|substitution_name, mutated_sites|
    site_list = SubstitutionEffects::SiteList.new(mutated_sites,
                                                fold_change_cutoff: fold_change_cutoff, pvalue_cutoff: pvalue_cutoff,
                                                to_be_disrupted: to_be_disrupted, to_be_preserved: to_be_preserved)
    [substitution_name, site_list]
  }.sort_by{|substitution_name, site_list|
    - site_list.quality
  }.reject{|substitution_name, site_list|
    output_configuration.suppress?(site_list)
  }.each{|substitution_name, site_list|
    seq_w_snv = sequence_with_snv_by_substitution_name(sequence, substitution_name)

    stream.puts ">#{substitution_name}"
    stream.puts seq_w_snv
    output_configuration.output(site_list, stream: stream)
    stream.puts
  }
end

fold_change_cutoff = 5.0
pvalue_cutoff = 0.0005

only_sites_overlapping_snv = false

output_configuration = SubstitutionEffects::OutputConfiguration.default

output_configuration.suppress_on_disrupted_sites_missing = false
output_configuration.suppress_on_preserved_sites_missing = false

output_configuration.suppress_on_disrupted_what_should_be_preserved = false
output_configuration.suppress_on_preserved_what_should_be_disrupted = false
output_configuration.suppress_on_emerged_what_should_be_disrupted = false

to_be_disrupted_patterns = []
to_be_preserved_patterns = []
OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <SNV sequence> [options]"
  opts.separator 'Options:'
  opts.on('--fold-change-cutoff CUTOFF', 'Fold change threshold to treat P-value change as significant') {|value|
    fold_change_cutoff = Float(value)
  }
  opts.on('--pvalue-cutoff CUTOFF', 'P-value to treat a word as a site') {|value|
    pvalue_cutoff = Float(value)
  }
  opts.on('--disrupt MOTIF_PATTERN', '-d', 'Specify motif name pattern for motifs to be disrupted. Can be used several times') {|value|
    to_be_disrupted_patterns << value
  }
  opts.on('--preserve MOTIF_PATTERN', '-p', 'Specify motif name pattern for motifs to be preserved. Can be used several times') {|value|
    to_be_preserved_patterns << value
  }

  opts.on('--only-sites-overlapping-snv') {|value|
    only_sites_overlapping_snv = true
  }

  opts.on('--hide-attentions', 'Don\'t show attentions about missing sites') {
    output_configuration.show_attentions = false
  }
  opts.on('--hide-strong-violations', 'Don\'t show strong violations about sites whose behavior is opposite to desired') {
    output_configuration.show_strong_violations = false
  }
  opts.on('--show-all-sites', 'Show all sites in sequence (not only sites of interest)') {
    output_configuration.show_all_sites = true
  }
  opts.on('--show-site-details', 'Show site details') {
    output_configuration.show_site_details = true
  }

  opts.on('--suppress-on-sites-missing', 'Suppress position when any site from targets (for disruption or preservation) is missing') {
    output_configuration.suppress_on_disrupted_sites_missing = true
    output_configuration.suppress_on_preserved_sites_missing = true
  }
  opts.on('--suppress-on-disrupted-sites-missing', 'Suppress position when any site from disruption targets is missing') {
    output_configuration.suppress_on_disrupted_sites_missing = true
  }
  opts.on('--suppress-on-preserved-sites-missing', 'Suppress position when any site from preservation targets is missing') {
    output_configuration.suppress_on_preserved_sites_missing = true
  }

  opts.on('--suppress-on-strong-violations', 'Suppress position when some criterium was strongly violated') {
    output_configuration.suppress_on_disrupted_what_should_be_preserved = true
    output_configuration.suppress_on_preserved_what_should_be_disrupted = true
    output_configuration.suppress_on_emerged_what_should_be_disrupted = true
  }

  opts.on('--suppress-on-disrupted-what-should-be-preserved', 'Suppress position when site which should be preserved was disrupted (strong violation)') {
    output_configuration.suppress_on_disrupted_what_should_be_preserved = true
  }

  opts.on('--suppress-on-preserved-what-should-be-disrupted', 'Suppress position when site which should be disrupted was preserved (strong violation)') {
    output_configuration.suppress_on_preserved_what_should_be_disrupted = true
  }

  opts.on('--suppress-on-emerged-what-should-be-disrupted', 'Suppress position when site which should be disrupted was emerged (strong violation)') {
    output_configuration.suppress_on_emerged_what_should_be_disrupted = true
  }

end.parse!(ARGV)

raise 'Specify SNV sequence'  unless sequence_with_snv = ARGV[0] # 'AAGCAGCGGCTTCTGAAGGAGGTAT[C/T]TATTTTGGTCCCAAACAGAAAAGAG'

sequence_with_snv = SequenceWithSNV.from_string(sequence_with_snv.upcase)

motif_names = Dir.glob('./motif_collection/*.pwm').map{|fn| File.basename(fn, '.pwm').to_sym }

to_be_disrupted = choose_motif_names_by_pattern(motif_names, to_be_disrupted_patterns)
to_be_preserved = choose_motif_names_by_pattern(motif_names, to_be_preserved_patterns)

$stderr.puts "Motifs to be disrupted:\n#{to_be_disrupted.to_a.join(',')}"
$stderr.puts "Motifs to be preserved:\n#{to_be_preserved.to_a.join(',')}"


sequence_with_snv.allele_variants.each_index do |i|
  sequence = sequence_with_snv.sequence_variant(i)
  sites_by_substitution = mutated_sites_by_substitution_name(sequence, pvalue_cutoff: pvalue_cutoff, ignore_positions: [sequence_with_snv.snv_position])

  if only_sites_overlapping_snv
    sites_by_substitution = select_sites_ovelapping_snv(sites_by_substitution, sequence_with_snv.left.length)
  end

  File.open("output_allele_#{sequence_with_snv.allele_variants[i]}.txt", 'w') do |fw|
    output_summary_of_sites(sequence, sites_by_substitution,
                            fold_change_cutoff: fold_change_cutoff, pvalue_cutoff: pvalue_cutoff,
                            to_be_disrupted: to_be_disrupted, to_be_preserved: to_be_preserved,
                            output_configuration: output_configuration,
                            stream: fw)
  end
end
