const TEMP_CACHE_FILE_PREFIX = "temp_cache"
CACHE_DIR = pwd()
DATA_KEY = :dat

function set_cache_dir(cache_dir)
    !isdir(cache_dir) && error(cache_dir, " not found!!!")
    global CACHE_DIR = cache_dir
end

temp_cache_file(state, cache_dir = CACHE_DIR) = 
    joinpath(cache_dir, "$(TEMP_CACHE_FILE_PREFIX)___$(hash(state)).bson")


function save_cache(state, data; cache_dir = CACHE_DIR,
        verbose = true, onerr::Function = (err) -> nothing)
    tcache_file = temp_cache_file(state, cache_dir) |> relpath
    try
        FileIO.save(tcache_file, Dict(DATA_KEY => data))
        verbose && tagprintln_inmw("CACHE SAVED\n", 
                "\ncache_file: ", tcache_file,
                "\nsize: ", filesize(tcache_file), " bytes",
                "\ndata type: ", typeof(data)
            )
    catch err
        verbose && tagprintln_inmw("ERROR SAVING CACHE\n", 
                "\ncache_file: ", tcache_file, 
                "\n", string_err(err)
            )

        onerr(err)
    end
end    

function load_cache(state; cache_dir = CACHE_DIR,
        verbose = true, onerr::Function = (err) -> nothing)

        tcache_file = temp_cache_file(state, cache_dir) |> relpath
    !isfile(tcache_file) && return nothing
    
    data = nothing
    try
        data = FileIO.load(tcache_file)[DATA_KEY]
        verbose && tagprintln_inmw("CACHE LOADED\n", 
                "\ncache_file: ", tcache_file,
                "\nsize: ", filesize(tcache_file), " bytes",
                "\ndata type: ", typeof(data)
            )
    catch err
        verbose && tagprintln_inmw("ERROR LOADING CACHE\n", 
                "\ncache_file: ", tcache_file, 
                "\n", string_err(err)
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
