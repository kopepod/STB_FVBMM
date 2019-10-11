classdef CODE 
 methods(Static = true)
 
  function X = codeFROMfile(dbdet,N,des,prefix)

   scr = '/local/java/ffmpeg/ffmpeg -i %s -vcodec rawvideo -y -loglevel quiet -pix_fmt rgb24 0.avi';
   
   X = cell(1e4,3);
   n = 0;
   
   c = cell(1,2);
   for k = 86
    load(strcat(dbdet,num2str(k),'.mat'));
    for j = 1:numel(c)
     if isempty(c{j}{2})
      continue
     end
     clc;disp(c{j}{1})
     v = VIDEO.readVideo(c{j}{1},scr,'0.avi');
     v = des.videoGsrc(v);
     c{j}{2} = CODE.codeVid(v,c{j}{2},des);
     if isempty(c{j}{2})
      continue
     end
     n = n+1;
     X{n,1} = c{j}{1};
     X{n,2} = c{j}{2};
     X{n,3} = k;
    end
    save(strcat(prefix,num2str(k),'.mat'),'c','-v7.3')
   end
   X = X(1:n,:);
  end
   
  
  function code = codeVid(v,stip,des)
   code = cell(size(stip,1),1);
   kpt  = true(size(stip,1),1);
   f_h = @des.encode;   
   parfor k = 1:size(stip,1)
    try
     code{k} = f_h(v,stip(k,:));
    catch
     kpt(k) = false;
    end
   end
   code = cell2mat(code(kpt));
   % code(end+1,1:3) = size(v);
   %figure
   %imagesc([code(:,4:end)])
  end
  
 end
end
