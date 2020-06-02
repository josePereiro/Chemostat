function search(model, hint; maxprint = 50, 
    fields = [:rxns, :mets, :rxnNames, :metNames, :genes])

hint = string(hint)
hint == "" && (println("0 found!!!"), return)

eqs_ = [] # equals
stw_ = [] # starts with
ctn_ = [] # contains
edw_ = [] # ends with

up = uppercase
up_hint = up(hint)
for field in fields
    dat = getfield(model, field)
    push!(eqs_, filter(x-> up(x) == up_hint, dat)...)
    push!(stw_, filter(x-> startswith(up(x), up_hint), dat)...)
    push!(ctn_, filter(x-> occursin(up_hint, up(x)), dat)...)
    push!(edw_, filter(x-> endswith(up(x), up_hint), dat)...)   
end

# print
all_res = [sort!(eqs_); sort!(stw_); sort!(ctn_); sort!(edw_)] |> unique
c = 0
println("$(length(all_res)) found!!!")
for res in all_res
    
    # I know maybe this is not the best way :)
    # This is only printing in bold the matches
    ci = 1 # current index
    while true
        hr_ = findnext(up_hint, up(res), ci)
        
        # No more hint
        if isnothing(hr_)
            print(res[ci:end])
            break;
        end
        
        # no-hint before hint
        print(res[ci:(first(hr_) - 1)])

        # hint
        printstyled(res[hr_], bold = true)

        ci = first(hr_) + length(hint)
        ci > length(res) && break
    end
    println()
    c == maxprint && (println("..."), return)
    c += 1
end
end