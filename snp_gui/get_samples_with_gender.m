function samples=get_samples_with_gender(samples, gender)
for i=1:length(samples)
    samples{i}=[samples{i} ' (' gender{i} ')'];
end