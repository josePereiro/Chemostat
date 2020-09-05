const TEMP_CACHE_FILE_PREFIX = "temp_cache"
CACHE_DIR = pwd()
DATA_KEY = :dat

function set_cache_dir(cache_dir)
    !isdir(cache_dir) && error(cache_dir, " not found!!!")
    global CACHE_DIR = cache_dir
end

temp_cache_file(state, cache_dir = CACHE_DIR, ext = ".jld") = 
    joinpath(cache_dir, string(TEMP_CACHE_FILE_PREFIX, "___", hash(state), ext))


function save_cache(state, data; cache_dir = CACHE_DIR, headline = "CACHE SAVED",
        verbose = true, onerr::Function = (err) -> rethrow(err))
    tcache_file = temp_cache_file(state, cache_dir) |> relpath
    try
        serialize(tcache_file, Dict(DATA_KEY => data))
        verbose && tagprintln_inmw(headline, 
                "\ncache_file: ", tcache_file,
                "\nsize: ", filesize(tcache_file), " bytes",
                "\ndata type: ", typeof(data),
                "\n"
            )
    catch err
        verbose && tagprintln_inmw("ERROR SAVING CACHE\n", 
                "\ncache_file: ", tcache_file, 
                "\n", string_err(err),
                "\n"
            )

        onerr(err)
    end
end    

function load_cache(state; cache_dir = CACHE_DIR, headline = "CACHE LOADED",
        verbose = true, onerr::Function = (err) -> rethrow(err))

        tcache_file = temp_cache_file(state, cache_dir) |> relpath
    !isfile(tcache_file) && return nothing
    
    data = nothing
    try
        data = deserialize(tcache_file)[DATA_KEY]
        verbose && tagprintln_inmw(headline, 
                "\ncache_file: ", tcache_file,
                "\nsize: ", filesize(tcache_file), " bytes",
                "\ndata type: ", typeof(data), 
                "\n"
            )
    catch err
        verbose && tagprintln_inmw("ERROR LOADING CACHE\n", 
                "\ncache_file: ", tcache_file, 
                "\n", string_err(err),
                "\n"
            )

        onerr(err)
    end
    return data
end    

function delete_temp_caches(cache_dir = CACHE_DIR; verbose = true)
    tcaches = filter(file -> startswith(file, TEMP_CACHE_FILE_PREFIX), readdir(cache_dir))
    for tc in tcaches
        tc = joinpath(cache_dir, tc)
        rm(tc, force = true)
        verbose && println(relpath(tc), " deleted!!!")
    end
end
