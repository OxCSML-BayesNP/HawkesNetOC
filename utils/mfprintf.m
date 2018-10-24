function mfprintf( fids, format, varargin )

% mfprintf writes tabular data to text file
% specified by fid
% in the required format 
% and with possible headers (specified in varargin)

for i=1:numel(fids)
    fprintf( fids(i), format, varargin{:} );
end

end