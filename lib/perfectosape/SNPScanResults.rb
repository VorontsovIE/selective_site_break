require 'set'

module PerfectosAPE
  Result = Struct.new(:line,
                      :variant_id, :motif_name,
                      :fold_change, :pvalue_1, :pvalue_2,
                      :pos_1, :orientation_1, :seq_1,
                      :pos_2, :orientation_2, :seq_2,
                      :variants ) do

    # "27610826_3 MAZ_f1  -7  direct  cggctgaGgaggaggag -7  direct  cggctgaCgaggaggag G/C 1.1218764110455249E-4 9.602413003842941E-4  0.11683275970285215"
    def self.from_string(line)
      variant_id, motif_name,
                pos_1, orientation_1, seq_1,
                pos_2, orientation_2, seq_2,
                variants,
                pvalue_1, pvalue_2, fold_change = line.split("\t")
      self.new( line,
                variant_id, motif_name.to_sym,
                fold_change.to_f, pvalue_1.to_f, pvalue_2.to_f,
                pos_1.to_i, orientation_1.to_sym, seq_1,
                pos_2.to_i, orientation_2.to_sym, seq_2,
                variants.to_sym )
    end

    def to_s
      line
    end

    def log2_fold_change
      Math.log2(fold_change)
    end

    def length
      seq_1.length
    end

    def allele_1
      letter = seq_1.each_char.detect{|letter| letter.upcase == letter }
      orientation_1 == :direct ? letter : Sequence.revcomp(letter)
    end

    def allele_2
      letter = seq_2.each_char.detect{|letter| letter.upcase == letter }
      orientation_2 == :direct ? letter : Sequence.revcomp(letter)
    end

    def best_allele
      fold_change < 1 ? allele_1 : allele_2
    end

    # Not used
    def seq_1_direct_strand
      orientation_1 == :direct ? seq_1 : Sequence.revcomp(seq_1)
    end

    def seq_2_direct_strand
      orientation_2 == :direct ? seq_2 : Sequence.revcomp(seq_2)
    end

    def seq_1_five_flank_length
      -pos_1
    end

    def seq_1_three_flank_length
      length - 1 + pos_1
    end

    def site_before_substitution?(pvalue_cutoff: 0.0005)
      pvalue_1 <= pvalue_cutoff
    end

    def site_after_substitution?(pvalue_cutoff: 0.0005)
      pvalue_2 <= pvalue_cutoff
    end

    def disrupted?(fold_change_cutoff: 5)
      fold_change <= (1.0 / fold_change_cutoff)
    end

    def emerged?(fold_change_cutoff: 5)
      fold_change >= fold_change_cutoff
    end

    def site_position_changed?
      pos_1 != pos_2 || orientation_1 != orientation_2
    end

    # glues corresponding motif site positions on direct and reverse strands
    def snv_position_in_site_1_pwm
      if orientation_1 == :direct
        - pos_1
      else
        pos_1 + length - 1
      end
    end

    def substitution_in_core?
      pos = snv_position_in_site_1_pwm
      pos >= 0 && pos < length
    end

    def substitution_in_flank?
      ! substitution_in_core?
    end

    # best site overlaps specified position
    def overlap_position?(pos)
      best_site_range.include?(pos)
    end

    def best_site_position
      (fold_change < 1) ? pos_1 : pos_2
    end

    def best_site_range
      best_site_position ... (best_site_position + length)
    end

    def best_site_word
      (fold_change < 1) ? seq_1 : seq_2
    end

    def worse_site_word
      (fold_change > 1) ? seq_1 : seq_2
    end

    # has site of specified strength at least on one allele
    def has_site_on_any_allele?(pvalue_cutoff:)
      site_before_substitution?(pvalue_cutoff: pvalue_cutoff) || site_after_substitution?(pvalue_cutoff: pvalue_cutoff)
    end

    def effect_strength_string(tabulated: false)
      if tabulated
        ['%7.2g' % log2_fold_change, '%7.2g' % pvalue_1, '%7.2g' % pvalue_2].join("\t")
      else
        "log2-Fold change %<fold_change>7.2g: from P-value of %<pvalue_1>7.2g to P-value of %<pvalue_2>7.2g" % {fold_change: log2_fold_change, pvalue_1: pvalue_1, pvalue_2: pvalue_2}
      end
    end

    def self.each_in_stream(stream, &block)
      stream.each_line.lazy.reject{|line|
        line.start_with?('#')
      }.map{|line|
        self.from_string(line)
      }.each(&block)
    end

    def self.each_in_file(all_mutations_filename, &block)
      return enum_for(:each_in_file, all_mutations_filename).lazy  unless block_given?
      File.open(all_mutations_filename) do |f|
        each_in_stream(f, &block)
      end
    end
  end
end
