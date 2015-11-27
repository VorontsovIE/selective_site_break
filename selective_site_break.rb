# Показать, какие P-value у каждого аллеля и fold-change
# Поставить галочки для мотивов нужного семейства и наоборот для вредных
# Показать сам сайт для каждого предсказания каждого фактора
# Писать отдельную секцию сайд-эффекты
# Ставить звездочку, если убит тот же самый сайт.
# показывать, что происходит с исходным мотивом
# показывать подпороги
# минимизировать число задетых семейств и показывать места, которые задевают мало других семейств (в случаях, где ничего не нашлось)

$:.unshift File.absolute_path('./lib', __dir__)
require 'perfectosape/SNPScanRunner'
require 'perfectosape/SNPScanResults'
require 'perfectosape/Hocomoco10Results'
require 'sequence_with_snv'
require 'site_list'
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


# SNV to string, marking reference position with lowercase letter
def format_snv(resulting_snv, pos_of_reference_snv)
  result = resulting_snv.to_s.upcase
  if pos_of_reference_snv < resulting_snv.left.length
    result[pos_of_reference_snv] = result[pos_of_reference_snv].downcase
  elsif pos_of_reference_snv == resulting_snv.left.length          
    result[pos_of_reference_snv, 5] = result[pos_of_reference_snv, 5].downcase
  else
    result[pos_of_reference_snv + 4] = result[pos_of_reference_snv + 4].downcase
  end
  result
end

def families_for_a_list_of_sites(sites)
  undefined_family_instead_of_empty = ->(families){ families.empty? ? ['Undefined'] : families }
  sites.map(&:motif_families).map(&undefined_family_instead_of_empty).flatten.uniq
end

############################################
def assess_mutations(site_infos, motif_families_to_disrupt,
                    position_to_overlap:,
                    fold_change_cutoff:, pvalue_cutoff:)
  disrupted_sites = site_infos.select{|site_info|
    site_info.disrupted?(fold_change_cutoff: fold_change_cutoff)
  }
  emerged_sites = site_infos.select{|site_info|
    site_info.emerged?(fold_change_cutoff: fold_change_cutoff)
  }
  # only those disrupted motifs which overlap specific position
  disrupted_sites_of_interest = disrupted_sites.select{|site_info|
    site_info.overlap_position?(position_to_overlap)
  }
  relocated_sites = site_infos.select{|site_info|
    site_info.has_site_on_any_allele?(pvalue_cutoff: pvalue_cutoff) && site_info.site_position_changed?
  }

  disrupted_motifs = disrupted_sites.map(&:motif_name).map(&:to_s)
  emerged_motifs = emerged_sites.map(&:motif_name).map(&:to_s)
  disrupted_motifs_of_interest = disrupted_sites_of_interest.map(&:motif_name).map(&:to_s)
  relocated_motifs = relocated_sites.map(&:motif_name).map(&:to_s)

  disrupted_motif_of_interest_families = disrupted_sites_of_interest.flat_map(&:motif_families).uniq
  disrupted_something_to_be_disrupted = disrupted_motif_of_interest_families.any?{|family|
    motif_families_to_disrupt.any?{|query_family| query_family == family }
  }
  affected_sites = (disrupted_sites + emerged_sites).uniq
  affected_families = families_for_a_list_of_sites(affected_sites)
  reliably_affected_families = families_for_a_list_of_sites(affected_sites.select(&:is_ABC_quality?))

  affect_something_not_to_be_affected = !(affected_families - motif_families_to_disrupt).empty?
  affect_something_reliable_not_to_be_affected = !(reliably_affected_families - motif_families_to_disrupt).empty?

  
  {
    # affected_motif_families: affected_motif_families,
    # affected_motifs: affected_motifs,
    disrupted_motifs: disrupted_motifs,
    emerged_motifs: emerged_motifs,
    relocated_motifs: relocated_motifs,
    disrupted_something_to_be_disrupted: disrupted_something_to_be_disrupted,
    affect_something_not_to_be_affected: affect_something_not_to_be_affected,
    affect_something_reliable_not_to_be_affected: affect_something_reliable_not_to_be_affected,
  }
end

def format_motif_name(motif_name)
  uniprot_id = motif_name.split('.').first
  families = WingenderTFClass::ProteinFamilyRecognizers::HumanAtLevel[3].subfamilies_by_uniprot_id(uniprot_id)
  "#{motif_name} (#{families.join(';')})"
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

motifs_to_disrupt_patterns = ARGV.drop(1).map{|pat| Regexp.new(pat, Regexp::IGNORECASE) }
motifs_to_disrupt = motif_names.select{|motif|
  motifs_to_disrupt_patterns.any?{|pat| pat.match(motif) }
}

motif_families_to_disrupt = motifs_to_disrupt.flat_map{|motif_name|
  uniprot_id = motif_name.to_s.split('.').first
  WingenderTFClass::ProteinFamilyRecognizers::HumanAtLevel[3].subfamilies_by_uniprot_id(uniprot_id)
}.uniq

puts "Motifs matching pattern: #{motifs_to_disrupt.join('; ')}"
puts "Motifs matching pattern: #{motif_families_to_disrupt.join('; ')}"
# motif_families_to_disrupt = ARGV.drop(1).map{|pat| Regexp.new(pat, Regexp::IGNORECASE)}
raise 'Specify motif families to disrupt'  if motif_families_to_disrupt.empty?
#############################################


# families_to_disrupt = WingenderTFClass::ProteinFamilyRecognizers::HumanAtLevel[3].subfamilies_by_uniprot_id(motif_uniprot_id_to_disrupt)

pos_of_reference_snv = sequence_with_snv.left.size

# # all_places of binding without check for P-value
# all_places = additionally_mutated_sites(sequence_with_snv)

sequence_with_snv.allele_variants.each_with_index{|allele, allele_index|
  puts "======================\nAllele variant: #{allele}\n======================"
  sequence = Sequence.new(sequence_with_snv.sequence_variant(allele_index))
  all_places = mutated_sites(sequence)

  # at least sites
  all_sites = all_places.select{|site_info|
    site_info.site_before_substitution?(pvalue_cutoff: pvalue_cutoff) ||
    site_info.site_after_substitution?(pvalue_cutoff: pvalue_cutoff)
  }

  all_sites.group_by(&:variant_id).each{|variant_id, site_infos|
    pos_of_snv = variant_id.sub(/^add:/, '').split(',').first.to_i
    result = assess_mutations(
      site_infos, motif_families_to_disrupt,
      position_to_overlap: pos_of_reference_snv - pos_of_snv,
      fold_change_cutoff: fold_change_cutoff, pvalue_cutoff: pvalue_cutoff
    )

    if !result[:affect_something_not_to_be_affected]
      status = "No sites of other TF families affected"
    elsif !result[:affect_something_reliable_not_to_be_affected]
      status = "Only bad-quality motifs among other TF families affected"
    else
      status = "Affect good-quality motifs of other TF families" # Doesn't go to output
    end

    if result[:disrupted_something_to_be_disrupted] && !result[:affect_something_reliable_not_to_be_affected]
      resulting_snv = sequence.sequence_with_snv_by_substitution_name(variant_id)
      resulting_snv_text = format_snv(resulting_snv, pos_of_reference_snv)

      puts "#{variant_id}\t#{status}\t#{resulting_snv_text}"
      unless result[:disrupted_motifs].empty? # Always true (for now)
        puts 'Disrupted: '
        puts result[:disrupted_motifs].map{|motif_name| "\t" + format_motif_name(motif_name) }.join("\n")
      end
      unless result[:emerged_motifs].empty?
        puts 'Emerged: '
        puts result[:emerged_motifs].map{|motif_name| "\t" + format_motif_name(motif_name) }.join("\n")
      end
      # $stderr.puts "Relocated: #{result[:relocated_motifs].join(' ')}"  unless result[:relocated_motifs].empty?
      # $stderr.puts '------------------------------------'
      puts
    end
  }
}
