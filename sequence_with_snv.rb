require_relative 'sequence'

class SequenceWithSNV
  PYRIMIDINES = ['C', 'T']
  attr_reader :left, :allele_variants, :right
  def initialize(left, allele_variants, right)
    # raise "SequenceWithSNV left part is invalid: #{left}" unless Sequence.valid_sequence?(left)
    # raise "SequenceWithSNV right part is invalid: #{right}" unless Sequence.valid_sequence?(right)
    # raise "SequenceWithSNV allele_variants are invalid: #{allele_variants}" unless allele_variants.map(&:to_s).all?{|letter| %w[A C G T N].include?(letter.upcase) }
    @left = left
    @allele_variants = allele_variants.map(&:to_s)
    @right = right
  end

  def self.from_string(sequence)
    left, mid, right = sequence.split(/[\[\]]/)
    allele_variants = mid.split('/')
    self.new(left, allele_variants, right)
  end

  def length
    left.length + 1 + right.length
  end

  def snv_position
    left.length
  end

  def revcomp
    SequenceWithSNV.new(Sequence.revcomp(right),
                        allele_variants.map{|letter| Sequence.complement(letter) },
                        Sequence.revcomp(left))
  end

  def in_pyrimidine_context?
    PYRIMIDINES.include?(allele_variants.first.upcase)
  end

  def in_pyrimidine_context
    in_pyrimidine_context? ? self : self.revcomp
  end

  def ==(other)
    other.is_a?(SequenceWithSNV) && @left == other.left && @allele_variants == other.allele_variants && @right == other.right
  end

  def eql?(other)
    other.class.equal?(SequenceWithSNV) && @left.eql?(other.left) && @allele_variants.eql?(other.allele_variants) && @right.eql?(other.right)
  end

  def hash
    [@left, @allele_variants, @right].hash
  end

  # main variant (or allele_variant_number variant) context
  def context(before: 1, after: 1, allele_variant_number: 0)
    left[-before..-1] + allele_variants[allele_variant_number] + right[0, after]
  end

  def subsequence(before:, after:)
    SequenceWithSNV.new(left[-before..-1], allele_variants, right[0, after])
  end

  def sequence_variant(allele_variant_number)
    "#{left}#{allele_variants[allele_variant_number]}#{right}"
  end

  def shuffle_string(str)
    str.each_char.to_a.shuffle.join
  end
  private :shuffle_string

  # shuffle flanks, but preserve 1bp-context
  def with_flanks_shuffled
    shuffled_left = shuffle_string(left[0..-2]) + left[-1] # preserve 1-bp context
    shuffled_right = right[0] + shuffle_string(right[1..-1])
    SequenceWithSNV.new(shuffled_left, allele_variants, shuffled_right)
  end

  def to_s
    allele_variants_str = allele_variants.join('/')
    "#{left}[#{allele_variants_str}]#{right}"
  end

  def inspect
    to_s
  end
end
