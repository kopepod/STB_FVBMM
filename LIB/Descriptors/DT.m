classdef DT
 methods(Static = true)
  
  function code = encode(file)
   
   ex = which('DenseTrack.');
   gz = '| gzip > tmp.gz';
   
   fileID = fopen('temp_script','w');
   fprintf(fileID,'%s %s %s',ex,which(file),gz);
   fclose(fileID);
   
   !echo -e '470' | sudo bash temp_script
   !gzip -d tmp.gz
   
   fileID = fopen('tmp','r');
   code = fscanf(fileID,'%f');
   fclose(fileID);
   
   code = code';
   
   code = reshape(code,40,[]);
   
   code = code';
   
   !rm tmp temp_script
   
  end
  
  function visualize(v,code)
   
   nv(:,:,1,:) = v;
   nv(:,:,2,:) = v;
   nv(:,:,3,:) = v;
   
   for k = 1:size(v,3)
    
    yx = round(code(code(:,1) == k,11:end));
    
    if ~isempty(yx)
     
     c = cell(size(yx,2)/2,1);
     m = 0;
     for j = 1:2:size(yx,2)
      m = m+1;
      c{m} = yx(:,j:j+1);
     end

     cl = rand(size(c{1},1),3);
     
     m = 0;
     for l = k-14:k
      m = m+1;

      Frame = nv(:,:,:,l);
      
      
      
      Frame = insertShape(Frame,...
       'Circle',[c{m},2*ones(size(c{m},1),1)],...
       'Color',cl,...
       'Opacity', 0.3);
      
      nv(:,:,:,l) = Frame;
      
     end
     
    end
   end
   
   implay(nv)
   
  end
  
 end
 
end