% function test_dna_utils()

AssignGeneralConstants;
%     AAGAAGACCCCAGAGAGAGACACACAGAGAGA                ATGACACACACGTGCT
%      AAAAAGACCCCAGAGAGAGACACACAGAGAGAAAAAAAAATAGTGAGTATGACACACACGTGCTGAC
s = ['AAGAAGACCCCAGAGAGAGACACACAGAGAGAGAGTTTTTGAGTGAGTATGACACACACGTGCTGAC'; ...
    'AAAAAAAAAAAAAGGGGGGGGGGGGATATGGGGGTGTGTGGGATGCGGGGGGAGAGTTTCCCCCCTC'];
%     AAAAAAAAAAAAAGGGGGGGGGGGGATATGGGAAAAATTGGGATGCGGGGGGAGAGTTTCCCCCCTC
%s = ['AAAAAAAACCCCCCCTTTTTTTT']; 

n= size(s,2);
[s_packed s_lens] = pack_seqs(s, matlab_word_size); 
s_again = int2nt(unpack_seqs(s_packed, s_lens, matlab_word_size))

for i=1:size(s,1)
   should_be_one = strcmp(s(i,:), s_again(i,:)) 
end



run_big=0;
if(run_big) % when we do run_big we still have a problem ..
    iters = 10000;
    s = ceil(rand(iters, n)*4);
else
    iters = 1;
end
[s_packed s_lens] = pack_seqs(s, matlab_word_size);
%s_packed = repmat(s_packed, iters, 1);
num_seqs = size(s,1);
L = 50; % This works now for L <= 16 and also L > 16
unique_flag = 1; % flag saying if we perform unique on the kmers
hash_flag = 0; % hash currently not supported  
int2nt(unpack_seqs(s_packed, s_lens, matlab_word_size));


kmers_packed = cell(2,1); kmer_inds = cell(2,1);
% s_double = zeros(size(s_packed), 'single');
% for i=1:size(s_packed,1)
%     s_double(i,:) = typecast(s_packed(i,:), 'single');
% end
% s_double = double(s_double);
for hash_flag = 0 % 1
    t_extract = cputime;
    [kmers_packed{hash_flag+1} kmer_inds{hash_flag+1}] = ...
        extract_sub_kmers(s_packed, repmat(s_lens, num_seqs, 1), L, unique_flag, hash_flag);
    
    [sorted_kmers_packed sort_perm] = sort(kmers_packed{hash_flag+1}); % sort kmers in increasing order
    sort_perm_inv = inv_perm(sort_perm);
    
    kmer_inds{hash_flag+1}(:,1) = sort_perm_inv(kmer_inds{hash_flag+1}(:,1));  % sort inds of kmers
    
    [aaa_val(hash_flag+1) aaa_ind(hash_flag+1)] = min(kmers_packed{hash_flag+1}(:,1));
    aaa_positions{hash_flag+1} = find(kmer_inds{hash_flag+1}(:,1) == aaa_ind(hash_flag+1));
    
    kmers = int2nt(unpack_seqs(kmers_packed{hash_flag+1}, L, matlab_word_size));
    S{hash_flag+1} =  sparse(kmer_inds{hash_flag+1}(:,1),kmer_inds{hash_flag+1}(:,2),1);
    %    figure; imagesc(S{hash_flag+1}); colorbar;
    t_extract = cputime - t_extract    
end
% figure; imagesc(S{2}-S{1}); colorbar;
if(length(S)>1)
    max_diff = max(abs(S{2}(:)-S{1}(:)))
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run first by input positions
kmer_seq_inds = [1 1 2 2 2 2 ]; 
kmer_positions_inds = [1 10 11 1 8 9]; 
num_kmers = length(kmer_seq_inds);
[kmers_packed_by_positions kmer_inds_by_positions] = ... % Another usage example: extract only some of the kmers: 
    extract_sub_kmers(s_packed, repmat(s_lens, num_seqs, 1), L, 0, hash_flag, ...
    kmer_seq_inds, kmer_positions_inds);


kmers_by_positions = int2nt(unpack_seqs(kmers_packed_by_positions, L, matlab_word_size))

for i=1:length(kmer_seq_inds)
    kmers_should_be(i,:) = s(kmer_seq_inds(i),kmer_positions_inds(i):kmer_positions_inds(i)+L-1);
end
kmers_should_be


ErrorStruct = []; ErrorStruct.error_model = 'exponential';
ErrorStruct.baseline_error = 0.01; 
ErrorStruct.final_error = 0.2; 
noise_table_vec = GenerateSubstitutionErrorTable(L, ErrorStruct);


seqs_weights_vec = GenerateMixtureWeights(num_seqs, 'power-law', 1); 
%seqs_weights_vec = ones(num_seqs,1)./num_seqs; % uniform frequency distribution 
num_reads = 1000; % number of reads to simulate 
[noisy_kmers_packed clean_kmers_packed kmer_inds] = ...
    SimulateReadsFromSequences(L, s_packed, s_lens, num_reads, ...
    seqs_weights_vec, noise_table_vec);

% % % noisy_kmers_packed = add_noise_to_kmers(L, num_kmers, ...
% % %     kmers_packed_by_positions, reshape(noise_table_vec, L, 16));
clean_kmers_by_positions = int2nt(unpack_seqs(clean_kmers_packed, L, matlab_word_size))
noisy_kmers_by_positions = int2nt(unpack_seqs(noisy_kmers_packed, L, matlab_word_size))
noise_matrix = num2str(noisy_kmers_by_positions ~= clean_kmers_by_positions) % matrix showing where simulated reads are different from original reads 


