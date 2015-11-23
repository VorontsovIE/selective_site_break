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


  # SNV constructors
  # sequence.make_seq_with_snv_by_position(10, 'A')
  def make_seq_with_snv_by_position(pos, substitution_nucleotide)
    raise 'Bad mutation nucleotide'  if substitution_nucleotide.to_s.upcase == sequence[pos].to_s.upcase
    SequenceWithSNV.new(sequence[0,pos], [sequence[pos], substitution_nucleotide], sequence[(pos + 1)..-1])
  end

  # sequence.sequence_with_snv_by_substitution_name('10,A')
  def sequence_with_snv_by_substitution_name(substitution_name)
    pos, substitution_nucleotide = substitution_name.chomp.split(',')
    make_seq_with_snv_by_position(pos.to_i, substitution_nucleotide)
  end

  def all_possible_snvs_at_position(pos)
    original_nucleotide = sequence[pos].upcase
    substitution_nucleotides = ['A','C','G','T'] - [original_nucleotide]
    substitution_nucleotides.map{|substitution_nucleotide|
      substitution_name = "#{pos},#{substitution_nucleotide}"
      seq_w_snv = make_seq_with_snv_by_position(pos, substitution_nucleotide)
      [substitution_name, seq_w_snv]
    }.to_h
  end

  def all_possible_snvs(ignore_positions: [])
    length.times.reject{|pos|
      ignore_positions.include?(pos)
    }.flat_map{|pos|
      all_possible_snvs_at_position(pos).to_a
    }.to_h
  end

end
