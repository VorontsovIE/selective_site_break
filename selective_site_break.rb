# (done) Показать, какие P-value у каждого аллеля и fold-change
# (done `+`/`-` и `*`) Поставить галочки для мотивов нужного семейства и наоборот для вредных
# (done) Показать сам сайт для каждого предсказания каждого фактора
# Писать отдельную секцию сайд-эффекты
# (done `!` и `?`) Ставить звездочку, если убит тот же самый сайт.
# (done) показывать, что происходит с исходным мотивом
# показывать подпороги
# минимизировать число задетых семейств и показывать места, которые задевают мало других семейств (в случаях, где ничего не нашлось)

$:.unshift File.absolute_path('./lib', __dir__)
require 'perfectosape/SNPScanRunner'
require 'perfectosape/SNPScanResults'
require 'perfectosape/Hocomoco10Results'
require 'sequence_with_snv'
require 'WingenderTFClass'

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

# Hocomoco10 intended
def all_sites_in_file(filename, additional_options: [])
  Ape.run_SNPScan(motif_collection: './motif_collection/',
                  snvs_file: filename,
                  precalulated_thresholds: './motif_thresholds',
                  fold_change_cutoff: 1,
                  pvalue_cutoff: 1,
                  additional_options: additional_options) do |pipe|
    PerfectosAPE::Hocomoco10Result.each_in_stream(pipe).to_a
  end
end

# snv_by_name is a hash {snv_name => seq_w_wnv}
def sites_for_several_snvs(snv_by_name)
  with_temp_file 'seq_generated.txt' do |f|
    snv_by_name.each{|substitution_name, seq_w_snv|
      f.puts "#{substitution_name}\t#{seq_w_snv}"
    }
    f.close
    all_sites_in_file(f.path)
  end
end

# Obtain affinity changes for every possible substitution in sequence
def mutated_sites(sequence, ignore_positions: [])
  sites_for_several_snvs(sequence.all_possible_snvs(ignore_positions: ignore_positions))
end

# Obtain affinity changes for SNV with additional nucleotide changed in a constant manner
def additionally_mutated_sites(sequence_with_snv)
  sites_for_several_snvs(sequence_with_snv.all_additional_substitutions)
end


def in_family?(site, families)
  !(site.motif_families & families).empty?
end

# SNV to string, marking reference position with lowercase letter
def format_snv(snv, pos_of_reference_snv)
  result = snv.to_s.upcase
  if pos_of_reference_snv < snv.left.length
    result[pos_of_reference_snv] = result[pos_of_reference_snv].downcase
  elsif pos_of_reference_snv == snv.left.length
    result[pos_of_reference_snv, 5] = result[pos_of_reference_snv, 5].downcase
  else
    result[pos_of_reference_snv + 4] = result[pos_of_reference_snv + 4].downcase
  end
  result
end

# TODO(refactor): move to a class
def families_for_a_list_of_sites(sites)
  undefined_family_instead_of_empty = ->(families){ families.empty? ? ['Undefined'] : families }
  sites.map(&:motif_families).map(&undefined_family_instead_of_empty).flatten.uniq
end

def families_by_motif_names(motif_names, level: 3)
  motif_names.flat_map{|motif_name|
    uniprot_id = motif_name.to_s.split('.').first
    WingenderTFClass::ProteinFamilyRecognizers::HumanAtLevel[level].subfamilies_by_uniprot_id(uniprot_id)
  }.uniq
end

############################################
class EffectAssessment
  attr_reader :sites, :fold_change_cutoff, :pvalue_cutoff
  # sites don't need to be actual sites on any allele
  def initialize(sites, fold_change_cutoff:, pvalue_cutoff:)
    @sites = sites
    @fold_change_cutoff = fold_change_cutoff
    @pvalue_cutoff = pvalue_cutoff
  end

  # sites having site on any allele
  def actual_sites
    @actual_sites ||= sites.select{|site|
      site.has_site_on_any_allele?(pvalue_cutoff: pvalue_cutoff)
    }
  end

  def disrupted_sites
    actual_sites.select{|site|
      site.disrupted?(fold_change_cutoff: fold_change_cutoff)
    }
  end

  def emerged_sites
    actual_sites.select{|site|
      site.emerged?(fold_change_cutoff: fold_change_cutoff)
    }
  end

  def relocated_sites
    actual_sites.select(&:site_position_changed?)
  end

  def affected_sites
    # (disrupted_sites + emerged_sites + relocated_sites).uniq
    (disrupted_sites + emerged_sites).uniq
  end

  def reliable_affected_sites
    affected_sites.select(&:is_ABC_quality?)
  end

  def affected_families
    families_for_a_list_of_sites(affected_sites)
  end

  def reliable_affected_families
    families_for_a_list_of_sites(reliable_affected_sites)
  end
end

# Effects for list of sites with regards to families of sites to be disrupted and the rest
class EffectAssessmentForSpecifiedFamilies < EffectAssessment
  attr_reader :motifs_to_disrupt, :position_to_overlap, :original_sites
  def initialize(sites, original_sites:, fold_change_cutoff:, pvalue_cutoff:, motifs_to_disrupt:, position_to_overlap:)
    super(sites, fold_change_cutoff: fold_change_cutoff, pvalue_cutoff: pvalue_cutoff)
    @motifs_to_disrupt = motifs_to_disrupt
    @position_to_overlap = position_to_overlap
    @original_sites = original_sites
  end

  def motif_families_to_disrupt
    @motif_families_to_disrupt ||= families_by_motif_names(motifs_to_disrupt)
  end

  # the same site which was broken by the reference SNP
  def original_site?(site)
    original_sites.select{|original_site|
      original_site.motif_name == site.motif_name
    }.select{|original_site|
      position_to_overlap + original_site.best_site_position == site.pos_1 ||
      position_to_overlap + original_site.best_site_position == site.pos_2
    }.size > 0
  end

  def disrupted_sites_of_interest
    disrupted_sites.select{|site_info|
      site_info.overlap_position?(position_to_overlap)
    }
  end

  def disrupted_something_to_be_disrupted?
    disrupted_motif_of_interest_families = disrupted_sites_of_interest.flat_map(&:motif_families).uniq
    disrupted_motif_of_interest_families.any?{|family|
      motif_families_to_disrupt.any?{|query_family| query_family == family }
    }
  end

  def affect_something_not_to_be_affected?
    !(affected_families - motif_families_to_disrupt).empty?
  end

  def affect_something_reliable_not_to_be_affected?
    !(reliable_affected_families - motif_families_to_disrupt).empty?
  end

  def status
    if !affect_something_not_to_be_affected?
      status = "No sites of other TF families affected"
    elsif !affect_something_reliable_not_to_be_affected?
      status = "Only bad-quality motifs among other TF families affected"
    else
      status = "Affect good-quality motifs of other TF families" # Doesn't go to output
    end
  end

  def site_list_formatted_string(list_of_sites, header:, indent: "\t")
    if list_of_sites.empty?
      nil
    else
      "#{header}:\n" + \
      list_of_sites.map{|site|
        infos = [
          site.motif_name_formatted,
          original_site?(site) ? '!' : '?',
          in_family?(site, motif_families_to_disrupt) ? '*' : '-',
          site.effect_strength_string,
          "#{site.seq_1} --> #{site.seq_2}",
        ]
        indent + infos.join("\t")  # + (motif_families_to_disrupt - site.motif_families).empty?
      }.join("\n") + "\n"
    end
  end
end

def process_snv(snv, variant_id, sites,
                stream: $stdout,
                original_sites:,
                motifs_to_disrupt:,
                pos_of_reference_snv:,
                fold_change_cutoff:, pvalue_cutoff:)
  pos_of_snv = variant_id.split(',').first.to_i
  result = EffectAssessmentForSpecifiedFamilies.new(sites,
    original_sites: original_sites,
    fold_change_cutoff: fold_change_cutoff, pvalue_cutoff: pvalue_cutoff,
    motifs_to_disrupt: motifs_to_disrupt,
    position_to_overlap: pos_of_reference_snv - pos_of_snv)

  if result.disrupted_something_to_be_disrupted? && !result.affect_something_reliable_not_to_be_affected?
    snv_text = format_snv(snv, pos_of_reference_snv)

    requested_to_disrupt = sites.select{|site| motifs_to_disrupt.include?(site.motif_name) }
    stream.puts "#{variant_id}\t#{result.status}\t#{snv_text}"
    stream.puts  result.site_list_formatted_string(requested_to_disrupt, header: "Requested to disrupt")
    stream.print result.site_list_formatted_string(result.disrupted_sites, header: "Disrupted")
    stream.print result.site_list_formatted_string(result.emerged_sites, header: "Emerged")
    stream.print result.site_list_formatted_string(result.relocated_sites, header: "Relocated")
    stream.puts
  end
end

#############################################
motif_names = Dir.glob('./motif_collection/*.pwm').map{|fn| File.basename(fn, '.pwm').to_sym }

fold_change_cutoff = 4.0
pvalue_cutoff = 0.0005
OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <SNV sequence> <patterns> [options]"
  opts.separator 'Options:'
  opts.on('--fold-change-cutoff CUTOFF', 'Fold change threshold to treat P-value change as significant') {|value|
    fold_change_cutoff = Float(value)
  }
  opts.on('--pvalue-cutoff CUTOFF', 'P-value to treat a word as a site') {|value|
    pvalue_cutoff = Float(value)
  }
end.parse!(ARGV)
# $stderr.puts(fold_change_cutoff: fold_change_cutoff, pvalue_cutoff: pvalue_cutoff)

raise 'Specify SNV sequence'  unless sequence_with_snv = ARGV[0] # 'AAGCAGCGGCTTCTGAAGGAGGTAT[C/T]TATTTTGGTCCCAAACAGAAAAGAG'
sequence_with_snv = SequenceWithSNV.from_string(sequence_with_snv.upcase)

patterns = ARGV.drop(1)
motifs_to_disrupt_patterns = patterns.map{|pat| Regexp.new(pat, Regexp::IGNORECASE) }
motifs_to_disrupt = motif_names.select{|motif|
  motifs_to_disrupt_patterns.any?{|pat| pat.match(motif) }
}

motif_families_to_disrupt = families_by_motif_names(motifs_to_disrupt)

puts "Motifs matching pattern: #{motifs_to_disrupt.join('; ')}"
puts "Motifs matching pattern: #{motif_families_to_disrupt.join('; ')}"
# motif_families_to_disrupt = ARGV.drop(1).map{|pat| Regexp.new(pat, Regexp::IGNORECASE)}
raise 'Specify motif families to disrupt'  if motif_families_to_disrupt.empty?
#############################################

puts sequence_with_snv
puts "Factor: #{patterns.join(', ')}"

def sites_sorted_by_relevance(site_list, motifs_to_disrupt, pvalue_cutoff:)
  motif_families_to_disrupt = families_by_motif_names(motifs_to_disrupt)
  target_sites = site_list.select{|site|
    in_family?(site, motif_families_to_disrupt)
  }
  target_sites.sort_by{|site| [site.pvalue_1, site.pvalue_2].min }
end

original_sites = sites_for_several_snvs({'Original-SNP' => sequence_with_snv})
target_sites_sorted = sites_sorted_by_relevance(original_sites, motifs_to_disrupt, pvalue_cutoff: pvalue_cutoff)


puts target_sites_sorted.map{|site|
  infos = [
          motifs_to_disrupt.include?(site.motif_name) ? '*' : '',
          site.has_site_on_any_allele?(pvalue_cutoff: pvalue_cutoff) ? '+' : '-',
          # in_family?(site, motif_families_to_disrupt) ? '*' : '-',
          site.motif_name_formatted,
          site.effect_strength_string,
          "#{site.seq_1} --> #{site.seq_2}",
        ]
  infos.join("\t")
}

sequence_with_snv.allele_variants.each_with_index{|allele, allele_index|
  puts "======================\nAllele variant: #{allele}\n======================"
  sequence = Sequence.new(sequence_with_snv.sequence_variant(allele_index))
  all_places = mutated_sites(sequence)
  all_sites = all_places #.select{|site_info|
#    site_info.has_site_on_any_allele?(pvalue_cutoff: pvalue_cutoff)
#  }
  all_sites.group_by(&:variant_id).each{|variant_id, site_infos|
    snv = sequence.sequence_with_snv_by_substitution_name(variant_id)
    process_snv(snv, variant_id, site_infos,
                original_sites: original_sites,
                motifs_to_disrupt: motifs_to_disrupt,
                pos_of_reference_snv: sequence_with_snv.left.size,
                fold_change_cutoff: fold_change_cutoff, pvalue_cutoff: pvalue_cutoff)
  }
}
