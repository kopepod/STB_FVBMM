classdef VIDEO
  
  methods(Static = true)
    
    function v = readVideoSupported(FILE)
      
       
      VR = vision.VideoFileReader(FILE,'ImageColorSpace','Intensity');
      v = cell(1e4,1);
      t = 0;
      while ~isDone(VR)
        t = t+1;
        v{t} = step(VR);
      end
      release(VR);
      
      v = v(1:t);
      
      v = cat(3,v{:});
      
    end
    
    function v = readVideo(varargin)
      
      file = varargin{1};
      
      if nargin == 1
        videoFReader = vision.VideoFileReader(which(file),'ImageColorSpace','Intensity');
        v = cell(1,1e+3);
        n = 0;
        while ~isDone(videoFReader)
          n = n+1;v(n) = {step(videoFReader)};
        end
        v = v(1:n);
        release(videoFReader);
        v = reshape(cell2mat(v),[size(v{1}),n]);
        return
      end
      
      if nargin ~= 3
        return
      end
      
      unix_script = varargin{2};
      fileID = fopen('temp_script','w');
      fprintf(fileID,unix_script,which(file));
      fclose(fileID);
      !bash temp_script
      
      videoFReader = vision.VideoFileReader('0.avi','ImageColorSpace','Intensity');
      
      v = cell(1,1e+3);
      n = 0;
      while ~isDone(videoFReader)
        n = n+1;v{n} = step(videoFReader);
      end
      v = v(1:n);
      release(videoFReader);
      v = reshape(cell2mat(v),[size(v{1}),n]);
      !rm temp_script 0.avi
    end
    
    function v = vid2of(v,alg,vis)
      % convert video to optical flow
      if vis
        figure
      end
      
      for t = 1:size(v,3)
        if vis
          imshow(uint8(255 * v(:,:,t)))
          hold on
        end
        flow = estimateFlow(alg,v(:,:,t));
        v(:,:,t) = flow.Vx + 1j * flow.Vy;
        if vis
          quiver(real(v(:,:,t)),imag(v(:,:,t)))
          pause(0.5)
        end
      end
    end
    
    
    function v = visTraj(v,T)
      
      
      nv = cell(size(v,3),1);
      
      for t = 1:size(v,3)
        I = v(:,:,t);
        nv{t} = cat(3,I,I,I);
      end
      
      v = cat(4,nv{:});
      
      
      for k = 1:size(T,1)
        
        t = T(k,:);
        
        t0 = t(1)-14;
        
        t = shiftdim(t(2:end),1);
        t = [t(1:2:end),t(2:2:end)];
        
        
        for j = 1:15
          v(:,:,:,t0) = insertMarker(v(:,:,:,t0),t(j,:),'+','color','green','size',2);
          t0 = t0+1;
        end
        
        
      end
      
    end
    
  end
end