function true_partition =relabel(data_p, n_obs)
     true_partition = zeros(n_obs,1);
     true_partition(1) = 1;
     count = 1;
     for i = 2:n_obs
         pos = findfirst(data_p(1:i), data_p(i));
         if pos == i
             count = count + 1;
             true_partition(i) = count;
         else
             true_partition(i) = true_partition(pos);
         end
     end
      true_partition;
end


       function pos = findfirst(vec, point)
       
       for i=1:length(vec)
            if  point==vec(i)
                pos=i;
                break
            end
       end
       
       end