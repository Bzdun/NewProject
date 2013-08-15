library EnctLib;

function Enct(i : Integer) : Integer; stdcall;
begin
  Enct := i + 10;
end;

exports
        Enct;
end.