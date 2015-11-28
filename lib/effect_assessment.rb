class EffectAssessment
  attr_reader :sites, :fold_change_cutoff, :pvalue_cutoff, :strong_pvalue_cutoff
  # sites don't need to be actual sites on any allele
  def initialize(sites, fold_change_cutoff:, pvalue_cutoff:, strong_pvalue_cutoff:)
    @sites = sites
    @fold_change_cutoff = fold_change_cutoff
    @pvalue_cutoff = pvalue_cutoff
    @strong_pvalue_cutoff = strong_pvalue_cutoff
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

  def families_for_a_list_of_sites(sites)
    undefined_family_instead_of_empty = ->(families){ families.empty? ? ['Undefined'] : families }
    sites.map(&:motif_families).map(&undefined_family_instead_of_empty).flatten.uniq
  end
end
