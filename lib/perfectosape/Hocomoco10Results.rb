require 'set'
require_relative 'SNPScanResults'

module PerfectosAPE
  class Hocomoco10Result < Result

    def uniprot_id
      motif_name.to_s.split('.').first
    end

    def motif_quality
      motif_name.to_s.split('.').last
    end

    def is_ABC_quality?
      ['A','B','C'].include?(motif_quality)
    end

    def motif_families(level: 3)
      WingenderTFClass::ProteinFamilyRecognizers::HumanAtLevel[level].subfamilies_by_uniprot_id(uniprot_id)
    end

    def motif_name_formatted(level: 3)
      families_str = motif_families(level: level).join(';')
      "#{motif_name} (#{families_str})"
    end

    def self.from_generic_result(result)
      self.new(result.line,
              result.variant_id, result.motif_name,
              result.fold_change, result.pvalue_1, result.pvalue_2,
              result.pos_1, result.orientation_1, result.seq_1,
              result.pos_2, result.orientation_2, result.seq_2,
              result.variants)
    end

    def self.from_string(line)
      res = super(line)
      self.from_generic_result(res)
    end

    def self.each_in_stream(stream, &block)
      PerfectosAPE::Result.each_in_stream(stream).map{|res| self.from_generic_result(res) }.each(&block)
    end

    def self.each_in_file(filename, &block)
      PerfectosAPE::Result.each_in_file(filename).map{|res| self.from_generic_result(res) }.each(&block)
    end
  end
end
