function string2var(s::String,v::Any)
    s=symbol(s)
    @eval (($s) = ($v))
end
