module SubstitutionEffects
  OutputConfiguration = Struct.new(
    :show_all_sites, :hide_non_sites, :hide_not_changed, :show_attentions, :show_strong_violations, :show_site_details,
    :suppress_on_disrupted_sites_missing, :suppress_on_preserved_sites_missing,
    :suppress_on_disrupted_what_should_be_preserved, :suppress_on_preserved_what_should_be_disrupted,
    :suppress_on_emerged_what_should_be_disrupted
  ) do
    def self.default
      self.from_kwargs(
        show_all_sites: false,
        hide_non_sites: false,
        hide_not_changed: false,
        show_attentions: true,
        show_strong_violations: true,
        show_site_details: false,

        suppress_on_disrupted_sites_missing: false,
        suppress_on_preserved_sites_missing: false,

        suppress_on_disrupted_what_should_be_preserved: false,
        suppress_on_preserved_what_should_be_disrupted: false,
        suppress_on_emerged_what_should_be_disrupted: false
      )
    end

    def self.from_kwargs(
      show_all_sites:, hide_non_sites:, hide_not_changed:, show_attentions:, show_strong_violations:, show_site_details:,
      suppress_on_disrupted_sites_missing:, suppress_on_preserved_sites_missing:,
      suppress_on_disrupted_what_should_be_preserved:, suppress_on_preserved_what_should_be_disrupted:,
      suppress_on_emerged_what_should_be_disrupted:)
      self.new(
        show_all_sites, hide_non_sites, hide_not_changed, show_attentions, show_strong_violations, show_site_details,
        suppress_on_disrupted_sites_missing, suppress_on_preserved_sites_missing,
        suppress_on_disrupted_what_should_be_preserved,  suppress_on_preserved_what_should_be_disrupted,
        suppress_on_emerged_what_should_be_disrupted
      )
    end

    def suppress?(site_list)
      return true  if suppress_on_disrupted_sites_missing && !site_list.missing_from_disrupted.empty?
      return true  if suppress_on_preserved_sites_missing && !site_list.missing_from_preserved.empty?

      return true  if suppress_on_disrupted_what_should_be_preserved && !site_list.should_be_disrupted_but_preserved.empty?
      return true  if suppress_on_preserved_what_should_be_disrupted && !site_list.should_be_preserved_but_disrupted.empty?
      return true  if suppress_on_emerged_what_should_be_disrupted && !site_list.should_be_disrupted_but_emerged.empty?

      false
    end


    def motifs_string(msg, motifs)
      "#{msg}: #{motifs.to_a.sort.join(',')}"
    end
    private :motifs_string

    def output_motifs(msg, motifs, stream: $stdout)
      stream.puts(motifs_string(msg, motifs))  unless motifs.empty?
    end
    private :output_motifs

    def output(site_list, stream: $stdout)
      stream.puts "Quality: #{site_list.quality}"
      # Sites of interest
      output_motifs 'Preserved sites of interest', site_list.preserved_motifs_of_interest, stream: stream
      output_motifs 'Disrupted sites of interest', site_list.disrupted_motifs_of_interest, stream: stream
      output_motifs 'Emerged sites of interest', site_list.emerged_motifs_of_interest, stream: stream

      # Sites of all motifs and all disruption/emergence events (not only motifs of interest)
      if show_all_sites
        output_motifs 'Sites before substitution (all)', site_list.all_sites_before_substitution, stream: stream
        output_motifs 'Sites after substitution (all)', site_list.all_sites_after_substitution, stream: stream
        output_motifs 'Disrupted sites (all)', site_list.disrupted_motifs, stream: stream
        output_motifs 'Emerged sites (all)', site_list.emerged_motifs, stream: stream
      end

      # Attentions! Site which should be disrupted/emerged do not exist at all
      if show_attentions
        output_motifs 'Attention! Site should be disrupted, but there\'s no site', site_list.missing_from_disrupted, stream: stream
        output_motifs 'Attention! Site should be preserved, but there\'s no site', site_list.missing_from_preserved, stream: stream
      end

      # Event works opposite to desired
      if show_strong_violations
        output_motifs 'Strong violations! Site should be preserved, but was disrupted', site_list.should_be_preserved_but_disrupted, stream: stream
        output_motifs 'Strong violations! Site should be disrupted, but was preserved', site_list.should_be_disrupted_but_preserved, stream: stream
        output_motifs 'Strong violations! Site should be disrupted, but emerged', site_list.should_be_disrupted_but_emerged, stream: stream
      end

      # sites of interest detailed info (PerfectosAPE output)
      if show_site_details
        interest_sites = show_all_sites ? site_list.mutated_sites : site_list.mutated_sites_of_interest

        interest_sites.each do |site|
          is_site_before = site.site_before_substitution?(pvalue_cutoff: site_list.pvalue_cutoff)
          is_site_after = site.site_after_substitution?(pvalue_cutoff: site_list.pvalue_cutoff)

          has_disrupted = site.disrupted?(fold_change_cutoff: site_list.fold_change_cutoff) && is_site_before
          has_emerged = site.emerged?(fold_change_cutoff: site_list.fold_change_cutoff) && is_site_after

          next  if hide_non_sites && !(is_site_before || is_site_after)
          next  if hide_not_changed && !(has_disrupted || has_emerged)

          substitution_pos = site.variant_id.split(",").first.to_i
          site_msg = ''
          if is_site_before
            site_msg << "site before"
          else
            site_msg << "no site before"
          end

          if is_site_after
            site_msg << "; site after"
          else
            site_msg << "; no site after"
          end          

          if has_disrupted
            site_msg << "; disrupted"
          elsif site.emerged?(fold_change_cutoff: site_list.fold_change_cutoff) && is_site_after
            site_msg << "; emerged"
          else
            site_msg << "; not changed significantly"
          end


          infos = [
            site.variant_id, site.motif_name, 
            '%.2g' % site.pvalue_1, '%.2g' % site.pvalue_2, '%.2g' % site.fold_change,
            site_msg,
            substitution_pos + site.pos_1, site.orientation_1, site.seq_1, 
            substitution_pos + site.pos_2, site.orientation_2, site.seq_2,
          ]
          stream.puts infos.join("\t")
        end
      end
    end
  end
end
