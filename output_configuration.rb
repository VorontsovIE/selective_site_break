module SubstitutionEffects
  OutputConfiguration = Struct.new(
    :show_all_sites, :show_attentions, :show_strong_violations, :show_site_details,
    :suppress_on_disrupted_sites_missing, :suppress_on_preserved_sites_missing,
    :suppress_on_disrupted_what_should_be_preserved, :suppress_on_preserved_what_should_be_disrupted,
    :suppress_on_emerged_what_should_be_disrupted
  ) do
    def self.default
      self.from_kwargs(
        show_all_sites: false,
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
      show_all_sites:, show_attentions:, show_strong_violations:, show_site_details:,
      suppress_on_disrupted_sites_missing:, suppress_on_preserved_sites_missing:,
      suppress_on_disrupted_what_should_be_preserved:, suppress_on_preserved_what_should_be_disrupted:,
      suppress_on_emerged_what_should_be_disrupted:)
      self.new(
        show_all_sites, show_attentions, show_strong_violations, show_site_details,
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
        output_motifs 'Sites before substitution', site_list.all_sites_before_substitution, stream: stream
        output_motifs 'Sites after substitution', site_list.all_sites_after_substitution, stream: stream
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
        stream.puts(interest_sites)
      end
    end
  end
end
