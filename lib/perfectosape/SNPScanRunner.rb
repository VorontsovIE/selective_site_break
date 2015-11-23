module Ape
  def self.run_SNPScan(snvs_file:,
                      motif_collection:,                      

                      precalulated_thresholds: nil,
                      mode: 'mono',
                      fold_change_cutoff: nil,
                      pvalue_cutoff: nil,
                      additional_options: [],
                      jvm_options: [], # ['-Xmx1G']
                      output_file: nil,
                      &block
                      )
    case mode
    when /^mono$/i
      package = 'ru.autosome.perfectosape.SNPScan'
    when /^di$/i
      package = 'ru.autosome.perfectosape.di.SNPScan'
    else
      raise "Unknown mode `#{mode}`"
    end

    cmd = ['java', *jvm_options, '-cp', 'ape.jar', package]

    args = [motif_collection, snvs_file]

    opts = []
    opts += ['--precalc', precalulated_thresholds]  if precalulated_thresholds
    opts += ['--fold-change-cutoff', fold_change_cutoff.to_s]  if fold_change_cutoff
    opts += ['--pvalue-cutoff', pvalue_cutoff.to_s]  if pvalue_cutoff
    opts += additional_options

    if output_file
      $stderr.puts 'Warning! In run_SNPScan both `output_file` and block provided. Block will be ignored.'  if block_given?
      system(*cmd, *args, *opts, {out: output_file})
    else
      IO.popen([*cmd, *args, *opts], &block)
    end

  end
end
