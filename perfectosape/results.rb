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

    def normalized_snv_name
      variant_id.split("_", 2)[0]
    end

    def length
      seq_1.length
    end

    # Not used
    def seq_1_direct_strand
      orientation_1 == :direct ? seq_1 : Sequence.revcomp(seq_1)
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
