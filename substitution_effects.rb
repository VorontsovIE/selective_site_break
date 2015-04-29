class SubstitutionEffects
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
    result -= 0.5 * ((missing_from_disrupted.size.to_f / to_be_disrupted.size) + (missing_from_preserved.size.to_f / to_be_preserved.size))  # ratios
    # result -= 0.5 * (should_be_preserved_but_disrupted.size + should_be_disrupted_but_preserved.size + should_be_disrupted_but_emerged.size) # counts
    result
  end


  def motifs_string(msg, motifs)
    "#{msg}: #{motifs.to_a.sort.join(',')}"
  end

  def output_motifs(stream, msg, motifs)
    stream.puts motifs_string(msg, motifs)  unless motifs.empty?
  end

  def output(stream, show_all_sites: true, show_attentions: true, show_strong_violations: true, show_site_details: false)
    stream.puts "Quality: #{quality}"
    # Sites of interest
    output_motifs stream, 'Preserved sites of interest', preserved_motifs & motifs_of_interest
    output_motifs stream, 'Disrupted sites of interest', disrupted_motifs & motifs_of_interest
    output_motifs stream, 'Emerged sites of interest', emerged_motifs & motifs_of_interest

    # Sites of all motifs and all disruption/emergence events (not only motifs of interest)
    if show_all_sites
      output_motifs stream, 'Sites before substitution', all_sites_before_substitution
      output_motifs stream, 'Sites after substitution', all_sites_after_substitution
      output_motifs stream, 'Disrupted sites (all)', disrupted_motifs
      output_motifs stream, 'Emerged sites (all)', emerged_motifs
    end

    # Attentions! Site which should be disrupted/emerged do not exist at all
    if show_attentions
      output_motifs stream, 'Attention! Site should be disrupted, but there\'s no site', missing_from_disrupted
      output_motifs stream, 'Attention! Site should be preserved, but there\'s no site', missing_from_preserved
    end

    # Event works opposite to desired
    if show_strong_violations
      output_motifs stream, 'Strong violations! Site should be preserved, but was disrupted', should_be_preserved_but_disrupted
      output_motifs stream, 'Strong violations! Site should be disrupted, but was preserved', should_be_disrupted_but_preserved
      output_motifs stream, 'Strong violations! Site should be disrupted, but emerged', should_be_disrupted_but_emerged
    end

    # sites of interest detailed info (PerfectosAPE output)
    if show_site_details
      interest_sites = show_all_sites ? mutated_sites : mutated_sites.select{|site| motifs_of_interest.include?(site.motif_name) }
      stream.puts(interest_sites)
    end
  end
end
