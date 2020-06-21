function search(model, hint, fun = (x) -> false; maxprint = 50, 
    fields = [:rxns, :mets, :rxnNames, :metNames, :genes])

hint = string(hint)
hint == "" && (println("0 found!!!"), return)

fun_ = [] # custom filter
eqs_ = [] # equals
stw_ = [] # starts with
ctn_ = [] # contains
edw_ = [] # ends with

up = uppercase
up_hint = up(hint)
for field in fields
    dat = getfield(model, field)
    push!(fun_, map(idx -> (dat[idx], idx), findall(fun, dat))...)
    push!(eqs_, map(idx -> (dat[idx], idx), findall(x-> up(x) == up_hint, dat))...)
    push!(stw_, map(idx -> (dat[idx], idx), findall(x-> startswith(up(x), up_hint), dat))...)
    push!(ctn_, map(idx -> (dat[idx], idx), findall(x-> occursin(up_hint, up(x)), dat))...)
    push!(edw_, map(idx -> (dat[idx], idx), findall(x-> endswith(up(x), up_hint), dat))...)
end

# print
all_res = [sort!(fun_); sort!(eqs_); sort!(stw_); sort!(ctn_); sort!(edw_)] |> unique
c = 0
println("$(length(all_res)) found!!!")
for (res, idx) in all_res
    
    print(idx, ": ")
    _print_bold(res, hint)
    
    println()
    c == maxprint && (println("..."), return)
    c += 1
end
end

function _print_bold(text, in_bold)
    up_in_bold = uppercase(in_bold)
    up_text = uppercase(text)

    # I know maybe this is not the best way :)
    # This is only printing in bold the matches
    ci = 1 # current index
    while true
        hr_ = findnext(up_in_bold, up_text, ci)
        
        # No more hint
        if isnothing(hr_)
            print(text[ci:end])
            break;
        end
        
        # no-hint before hint
        print(text[ci:(first(hr_) - 1)])

        # hint
        printstyled(text[hr_], bold = true)

        ci = first(hr_) + length(in_bold)
        ci > length(text) && break
    end
end