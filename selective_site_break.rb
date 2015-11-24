$:.unshift File.absolute_path('./lib', __dir__)
require 'perfectosape/SNPScanRunner'
require 'perfectosape/SNPScanResults'
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

def all_sites_in_file(filename, additional_options: [])
  Ape.run_SNPScan(motif_collection: './motif_collection/',
                  snvs_file: filename,
                  precalulated_thresholds: './motif_thresholds',
                  fold_change_cutoff: 1,
                  pvalue_cutoff: 1,
                  additional_options: additional_options) do |pipe|
    PerfectosAPE::Result.each_in_stream(pipe).to_a
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
$stderr.puts(fold_change_cutoff: fold_change_cutoff, pvalue_cutoff: pvalue_cutoff)

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

$stderr.puts "Motifs matching pattern: #{motifs_to_disrupt.join('; ')}"
$stderr.puts "Motifs matching pattern: #{motif_families_to_disrupt.join('; ')}"
# motif_families_to_disrupt = ARGV.drop(1).map{|pat| Regexp.new(pat, Regexp::IGNORECASE)}
raise 'Specify motif families to disrupt'  if motif_families_to_disrupt.empty?
#############################################


# families_to_disrupt = WingenderTFClass::ProteinFamilyRecognizers::HumanAtLevel[3].subfamilies_by_uniprot_id(motif_uniprot_id_to_disrupt)

# all_places of binding without check for P-value
all_places = additionally_mutated_sites(sequence_with_snv)

# at least sites
all_sites = all_places.select{|site_info|
  site_info.site_before_substitution?(pvalue_cutoff: pvalue_cutoff) ||
  site_info.site_after_substitution?(pvalue_cutoff: pvalue_cutoff)
}

all_sites.group_by(&:variant_id).each{|variant_id, site_infos|
  disrupted_motifs = site_infos.select{|site_info|
    site_info.disrupted?(fold_change_cutoff: fold_change_cutoff)
  }.map(&:motif_name).map(&:to_s)

  emerged_motifs = site_infos.select{|site_info|
    site_info.emerged?(fold_change_cutoff: fold_change_cutoff)
  }.map(&:motif_name).map(&:to_s)

  relocated_motifs = site_infos.select{|site_info|
    (
      site_info.site_before_substitution?(pvalue_cutoff: pvalue_cutoff) ||
      site_info.site_after_substitution?(pvalue_cutoff: pvalue_cutoff)
    ) && site_info.site_position_changed?
  }.map(&:motif_name).map(&:to_s)


  disrupted_uniprot_ids = disrupted_motifs.map{|motif_name| motif_name.split('.').first }
  disrupted_motif_families = disrupted_uniprot_ids.flat_map{|uniprot_id|
    WingenderTFClass::ProteinFamilyRecognizers::HumanAtLevel[3].subfamilies_by_uniprot_id(uniprot_id)
  }.uniq

  # affected_motifs = (disrupted_motifs + emerged_motifs + relocated_motifs).uniq
  affected_motifs = (disrupted_motifs + emerged_motifs).uniq

  affected_uniprot_ids = affected_motifs.select{|motif_name|
      quality = motif_name.split('.').last
      ['A','B','C'].include?(quality) # Don't mind if D or S quality motif was affected
    }
    .map{|motif_name| motif_name.split('.').first }
  affected_motif_families = affected_uniprot_ids.flat_map{|uniprot_id|
    WingenderTFClass::ProteinFamilyRecognizers::HumanAtLevel[3].subfamilies_by_uniprot_id(uniprot_id)
  }.uniq


  disrupted_something_to_be_disrupted = disrupted_motif_families.any?{|family|
    motif_families_to_disrupt.any?{|query_family| query_family == family }
  }
  affect_something_not_to_be_affected = affected_motif_families.any?{|family|
    motif_families_to_disrupt.none?{|query_family| query_family == family }
  }

  if disrupted_something_to_be_disrupted && !affect_something_not_to_be_affected
    puts '------------------------------------'
    puts variant_id
    puts "Affected families: #{affected_motif_families.join('; ')}"
    puts "Disrupted: #{disrupted_motifs.join(' ')}"  unless disrupted_motifs.empty?
    puts "Emerged: #{emerged_motifs.join(' ')}"  unless emerged_motifs.empty?
    puts "Relocated: #{relocated_motifs.join(' ')}"  unless relocated_motifs.empty?
  end
}
