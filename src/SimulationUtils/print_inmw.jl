using Dates
using Distributed

MASTERW = 1
set_MASTERW(w) = (global MASTERW = w)

print_inmw(ss...) = (remotecall_wait(Core.print, MASTERW, ss...); nothing)
println_inmw(ss...) = (remotecall_wait(Core.println, MASTERW, ss...); nothing)
print_ifmw(ss...) = myid() == MASTERW ? Core.print(ss...) : nothing
println_ifmw(ss...) = myid() == MASTERW ? Core.println(ss...) : nothing


function wtag(io::IO)
    ws = min(displaysize(io) |> last, 80) 
    return lpad(
        string(
            " Worker ", myid(), " (", getpid(), ")",
            " [", Time(now()), "] "
        ), 
        ws, "-"
    )
end

function tagprint_inmw(ss...; tag::String = wtag(stdout))
    !isempty(tag) && println_inmw(tag)
    print_inmw(ss...)
end

function tagprintln_inmw(ss...; tag::String = wtag(stdout))
    !isempty(tag) && println_inmw(tag)
    println_inmw(ss...)
end

function tagprint_ifmw(ss...; tag::String = wtag(stdout))
    !isempty(tag) && println_ifmw(tag)
    print_ifmw(ss...)
end

function tagprintln_ifmw(ss...; tag::String = wtag(stdout))
    !isempty(tag) && println_ifmw(tag)
    println_ifmw(ss...)
end


function string_err(err; max_len = 10000)
    s = sprint(showerror, err, catch_backtrace())
    return length(s) > max_len ? s[1:max_len] * "\n..." : s
end