%function [CopyMatA CopyMatB] = replace_AB_according_to_strand(CopyMatA, CopyMatB, snp_ids, chip_type)
function [CopyMatA CopyMatB] = replace_AB_according_to_strand(CopyMatA, CopyMatB, snp_ids, chip_type)

ret_strand = load_snp_ids_strand(snp_ids, chip_type);

minus_strand = strmatch('-', ret_strand);
[CopyMatA(minus_strand, :) CopyMatB(minus_strand, :)] = swap(CopyMatA(minus_strand, :), CopyMatB(minus_strand, :));
