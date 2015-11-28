require_relative 'effect_assessment'

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
          EffectAssessmentForSpecifiedFamilies.in_family?(site, motif_families_to_disrupt) ? '*' : '-',
          site.effect_strength_string,
          "#{site.seq_1} --> #{site.seq_2}",
        ]
        indent + infos.join("\t")  # + (motif_families_to_disrupt - site.motif_families).empty?
      }.join("\n") + "\n"
    end
  end

  def self.in_family?(site, families)
    !(site.motif_families & families).empty?
  end
end
