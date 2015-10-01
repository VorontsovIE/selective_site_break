module SubstitutionEffects
  class SiteList
    attr_reader :mutated_sites, :fold_change_cutoff, :pvalue_cutoff, :to_be_disrupted, :to_be_preserved
    def initialize(mutated_sites,
                  fold_change_cutoff: , pvalue_cutoff:,
                  to_be_disrupted:, to_be_preserved:)
      @mutated_sites = mutated_sites
      @fold_change_cutoff = fold_change_cutoff
      @pvalue_cutoff = pvalue_cutoff
      @to_be_disrupted = to_be_disrupted
      @to_be_preserved = to_be_preserved
    end
    
    def motifs_of_interest
      to_be_disrupted | to_be_preserved
    end

    def mutated_sites_of_interest
      mutated_sites.select{|site|
        motifs_of_interest.include?(site.motif_name)
      }
    end

    # selects motifs by corresponding sites
    def select_motifs(&block)
      mutated_sites.select{|site_info|
        block.call(site_info)
      }.map(&:motif_name).to_set
    end

    def all_sites_before_substitution
      select_motifs{|site_info|
        site_info.site_before_substitution?(pvalue_cutoff: pvalue_cutoff)
      }
    end

    def all_sites_after_substitution
      select_motifs{|site_info|
        site_info.site_after_substitution?(pvalue_cutoff: pvalue_cutoff)
      }
    end
    
    # fold change check in perfectosape is bidirectional, so check again
    def disrupted_motifs
      select_motifs{|site_info|
        site_info.disrupted?(fold_change_cutoff: fold_change_cutoff)
      } & all_sites_before_substitution
    end

    def preserved_motifs
      select_motifs{|site_info|
        !site_info.disrupted?(fold_change_cutoff: fold_change_cutoff)
      } & all_sites_before_substitution
    end

    def emerged_motifs
      select_motifs{|site_info|
        site_info.emerged?(fold_change_cutoff: fold_change_cutoff)
      } & all_sites_after_substitution
    end

    # Reduction to motifs of interest
    def preserved_motifs_of_interest
      preserved_motifs & motifs_of_interest
    end

    def disrupted_motifs_of_interest
      disrupted_motifs & motifs_of_interest
    end

    def emerged_motifs_of_interest
      emerged_motifs & motifs_of_interest
    end

    # sites should be disrupted/preserved but sequence don't even have motif's site
    def missing_from_disrupted
      to_be_disrupted - all_sites_before_substitution
    end

    def missing_from_preserved
      to_be_preserved - all_sites_before_substitution
    end

    # sites which are really present but were disrupted
    def actual_sites_to_be_disrupted
      to_be_disrupted & all_sites_before_substitution
    end

    def actual_sites_to_be_preserved
      to_be_preserved & all_sites_before_substitution
    end

    # Very bad cases
    def should_be_preserved_but_disrupted
      disrupted_motifs & to_be_preserved
    end

    def should_be_disrupted_but_preserved
      preserved_motifs & to_be_disrupted
    end

    def should_be_disrupted_but_emerged
      emerged_motifs & to_be_disrupted
    end

    # Quality metrics
    def quality_disruption
      result = (disrupted_motifs & actual_sites_to_be_disrupted).size.to_f / (actual_sites_to_be_disrupted).size
      result.nan? ? 0 : result
    end
    def quality_preservation
      result = (preserved_motifs & actual_sites_to_be_preserved).size.to_f / (actual_sites_to_be_preserved).size
      result.nan? ? 0 : result
    end
    def quality
      result = 0
      result += 1.0 * quality_disruption  # ratio
      result += 1.0 * quality_preservation  # ratio
      missing_from_disrupted_part = (to_be_disrupted.size != 0) ? (missing_from_disrupted.size.to_f / to_be_disrupted.size) : 0
      missing_from_preserved_part = (to_be_preserved.size != 0) ? (missing_from_preserved.size.to_f / to_be_preserved.size) : 0
      result -= 0.5 * (missing_from_disrupted_part + missing_from_preserved_part)  # ratios
      # result -= 0.5 * (should_be_preserved_but_disrupted.size + should_be_disrupted_but_preserved.size + should_be_disrupted_but_emerged.size) # counts
      result
    end
  end
end
