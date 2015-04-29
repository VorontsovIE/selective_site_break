class Sequence
  attr_reader :sequence
  def initialize(sequence)
    # raise "Wrong sequence `#{sequence}`"  unless Sequence.valid_sequence?(sequence)
    @sequence = sequence
  end

  def length
    sequence.length
  end

  def revcomp
    Sequence.new(Sequence.revcomp(sequence))
  end

  def ==(other)
    other.is_a?(Sequence) && @sequence == other.sequence
  end

  def eql?(other)
    other.class.equal?(Sequence) && @sequence.eql?(other.sequence)
  end

  def hash
    @sequence.hash
  end

  ACGT = 'acgtACGT'.freeze
  TGCA = 'tgcaTGCA'.freeze
  def self.complement(sequence)
    sequence.tr(ACGT, TGCA)
  end

  def self.revcomp(sequence)
    complement(sequence).reverse!
  end

  def self.valid_sequence?(sequence)
    sequence.match /\A[acgtn]+\z/i
  end
end
