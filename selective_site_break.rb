$:.unshift File.absolute_path('./lib', __dir__)
require 'perfectosape/SNPScanRunner'
require 'perfectosape/SNPScanResults'
require 'perfectosape/Hocomoco10Results'
require 'sequence_with_snv'
require 'WingenderTFClass'

require 'effect_assessment_for_family'

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

def families_by_motif_names(motif_names, level: 3)
  motif_names.flat_map{|motif_name|
    uniprot_id = motif_name.to_s.split('.').first
    WingenderTFClass::ProteinFamilyRecognizers::HumanAtLevel[level].subfamilies_by_uniprot_id(uniprot_id)
  }.uniq
end

############################################
def process_snv(snv, variant_id,
                effect_assessment,
                stream: $stdout,
                pos_of_reference_snv:)

  if ! effect_assessment.desired_effects.empty? && \
      families_by_motif_names( effect_assessment.side_effects.select(&:is_ABC_quality?).map(&:motif_name) ).uniq.size <= 2 #&& \
      effect_assessment.side_effects.select(&:is_ABC_quality?).uniq.size <= 5 #!effect_assessment.affect_something_reliable_not_to_be_affected?
    snv_text = format_snv(snv, pos_of_reference_snv)

    stream.puts "#{variant_id}\t#{effect_assessment.status}\t#{snv_text}"
    stream.print effect_assessment.site_list_formatted_string(effect_assessment.desired_effects, snv, header: "Desired effects")
    stream.print effect_assessment.site_list_formatted_string(effect_assessment.side_effects, snv, header: "Side effects")
    stream.print effect_assessment.site_list_formatted_string(effect_assessment.on_edge_effects, snv, header: "Edge effects")
    # stream.puts  effect_assessment.site_list_formatted_string(effect_assessment.sites_requested_to_disrupt, snv, header: "Requested to disrupt")
    # stream.print effect_assessment.site_list_formatted_string(effect_assessment.disrupted_sites, snv, header: "Disrupted")
    # stream.print effect_assessment.site_list_formatted_string(effect_assessment.emerged_sites, snv, header: "Emerged")
    # stream.print effect_assessment.site_list_formatted_string(effect_assessment.relocated_sites, snv, header: "Relocated")
    stream.puts
  end
end


def sites_sorted_by_relevance(site_list, motifs_to_disrupt, pvalue_cutoff:)
  motif_families_to_disrupt = families_by_motif_names(motifs_to_disrupt)
  target_sites = site_list.select{|site|
    EffectAssessmentForSpecifiedFamilies.in_family?(site, motif_families_to_disrupt)
  }
  target_sites.sort_by{|site| [site.pvalue_1, site.pvalue_2].min }
end

#############################################
motif_names = Dir.glob('./motif_collection/*.pwm').map{|fn| File.basename(fn, '.pwm').to_sym }

fold_change_cutoff = 4.0
pvalue_cutoff = 0.001
strong_pvalue_cutoff = 0.0005
OptionParser.new do |opts|
  opts.banner = "Usage: #{opts.program_name} <SNV sequence> <patterns> [options]"
  opts.separator 'Options:'
  opts.on('--fold-change-cutoff CUTOFF', 'Fold change threshold to treat P-value change as significant') {|value|
    fold_change_cutoff = Float(value)
  }
  opts.on('--pvalue-cutoff CUTOFF', 'P-value to treat a word as a weak site') {|value|
    pvalue_cutoff = Float(value)
  }
  opts.on('--strong-pvalue-cutoff CUTOFF', 'P-value to treat a word as a strong site') {|value|
    strong_pvalue_cutoff = Float(value)
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

# motif_families_to_disrupt = ARGV.drop(1).map{|pat| Regexp.new(pat, Regexp::IGNORECASE)}
raise 'Specify motif families to disrupt'  if motif_families_to_disrupt.empty?

puts 'Legend:'
puts '  Motif qualifiers:'
puts '    ! - target motif'
puts '    ~ - motif of target family'
puts '    # - motif of another family'
puts '    * - motif is reliable (A/B/C quality in Hocomoco10)'
puts '  Site qualifiers:'
puts "    S - strong site (P-value #{strong_pvalue_cutoff})"
puts "    w - weak site (P-value #{pvalue_cutoff})"
puts "    n - not a site"
puts "    r - best position of motif was relocated"
puts '--------------------------------------------'

puts "Original SNP: #{sequence_with_snv}"
puts "Target transcription factor: #{patterns.join(', ')}"
puts "Target motif (to induce affinity loss) !: #{motifs_to_disrupt.join('; ')}"
puts "Target family ~: #{motif_families_to_disrupt.join('; ')}"
puts '--------------------------------------------'
#############################################

original_sites = sites_for_several_snvs({'Original-SNP' => sequence_with_snv})
target_sites_sorted = sites_sorted_by_relevance(original_sites, motifs_to_disrupt, pvalue_cutoff: pvalue_cutoff)
best_alleles = motifs_to_disrupt.map{|motif_name| original_sites.detect{|site| site.motif_name == motif_name }.best_allele }.uniq
if best_alleles.size == 1
  best_allele = best_alleles[0]
  puts "Base allele (with stronger prediction): #{best_allele}"
elsif best_alleles.size > 1
  puts "Diagnosis: several target motifs have different affinity change direction for a given SNP. Please check your input data!"
else
  raise 'WTF'
end

headers = [
  '', '', 'Motif', 'log2-Fold change',
  "P-value #{sequence_with_snv.allele_variants[0]}", "P-value #{sequence_with_snv.allele_variants[1]}",
  "Binding site #{sequence_with_snv.allele_variants[0]}", "Binding site #{sequence_with_snv.allele_variants[1]}",
  "Allele with stronger prediction", 'TF family'
]

effect_assessment_original = EffectAssessment.new(original_sites, fold_change_cutoff: fold_change_cutoff, pvalue_cutoff: pvalue_cutoff, strong_pvalue_cutoff: strong_pvalue_cutoff)
puts effect_assessment_original.site_list_formatted_string(target_sites_sorted, sequence_with_snv, motifs_to_disrupt, motif_families_to_disrupt, header: "Sites of family requested to disrupt overlapping reference SNP")
puts '--------------------------------------------------'
puts 'Suggested substitutions'
puts '--------------------------------------------------'

allele_index = sequence_with_snv.allele_variants.index{|allele| allele == best_allele }
sequence = Sequence.new(sequence_with_snv.sequence_variant(allele_index))
all_places = mutated_sites(sequence)
all_sites = all_places.select{|site_info|
   site_info.has_site_on_any_allele?(pvalue_cutoff: pvalue_cutoff)
 }
all_sites.group_by(&:variant_id).map{|variant_id, sites|
  snv = sequence.sequence_with_snv_by_substitution_name(variant_id)
  pos_of_reference_snv = sequence_with_snv.left.size

  pos_of_snv = variant_id.split(',').first.to_i
  effect_assessment = EffectAssessmentForSpecifiedFamilies.new(sites,
    original_sites: original_sites,
    motifs_to_disrupt: motifs_to_disrupt,
    position_to_overlap: pos_of_reference_snv - pos_of_snv,
    fold_change_cutoff: fold_change_cutoff, pvalue_cutoff: pvalue_cutoff, strong_pvalue_cutoff: strong_pvalue_cutoff)

  {effect_assessment: effect_assessment, snv: snv, variant_id: variant_id, pos_of_reference_snv: pos_of_reference_snv}
}.sort_by{|infos|
  infos[:effect_assessment].reliable_erroneously_affected_families.size
}.each{|infos|
  process_snv(infos[:snv], infos[:variant_id], infos[:effect_assessment], pos_of_reference_snv: infos[:pos_of_reference_snv])
}
