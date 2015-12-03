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

  def site_strength_code(site)
    (site.has_site_on_any_allele?(pvalue_cutoff: strong_pvalue_cutoff) ? 'S' : (site.has_site_on_any_allele?(pvalue_cutoff: pvalue_cutoff) ? 'w' : 'n')) +
            (site.site_position_changed? ? 'r' : '')
  end

  def site_relevance_code(site, motifs_to_disrupt, motif_families_to_disrupt)
    (motifs_to_disrupt.include?(site.motif_name) ? '!' : (in_family?(site, motif_families_to_disrupt) ? '~' : '#')) +
            (site.is_ABC_quality? ? '*' : '')
  end

  def site_infos(site, motifs_to_disrupt, motif_families_to_disrupt)
    [
      site_relevance_code(site, motifs_to_disrupt, motif_families_to_disrupt),
      site_strength_code(site),
      site.motif_name,
      site.effect_strength_string(tabulated: true),
      site.seq_1, site.seq_2,
      site.best_allele,
      site.motif_families.join('; '),
    ]
  end

  def site_list_formatted_string(list_of_sites, snv, motifs_to_disrupt, motif_families_to_disrupt, header:, indent: "")
    if list_of_sites.empty?
      nil
    else
      headers = [
        '', '', 'Motif', 'log2-Fold change',
        "P-value #{snv.allele_variants[0]}", "P-value #{snv.allele_variants[1]}",
        "Binding site #{snv.allele_variants[0]}", "Binding site #{snv.allele_variants[1]}",
        "Allele with stronger prediction", 'TF family'
      ]

      "#{header}:\n" + headers.join("\t") + "\n" + list_of_sites.map{|site|
        indent + site_infos(site, motifs_to_disrupt, motif_families_to_disrupt).join("\t")
      }.join("\n") + "\n"
    end
  end

  def in_family?(site, families)
    EffectAssessment.in_family?(site, families)
  end

  def self.in_family?(site, families)
    !(site.motif_families & families).empty?
  end
end
