# The color plan below is better in terminal dark mode
"""
Print a title.
"""
function title(msg)
    !echo_on && return
    printstyled('\n', msg, '\n', bold = true, color = :cyan)
    printstyled(repeat('=', length(msg) + 2), '\n', color=31)
end

"""
Print a subtitle.
"""
function subtitle(msg)
    !echo_on && return
    printstyled('\n', msg, '\n', color = :cyan)
    printstyled(repeat('-', length(msg) + 2), '\n', color=31)
end

"""
Print a message.
"""
function message(msg...)
    !echo_on && return
    printstyled("\nMessage: \n", color = :light_magenta)
    for m in msg
        printstyled("    - $m\n", color = 86)
    end
end

"""
Print a warning.
"""
function warning(msg...)
    !echo_on && return
    printstyled("\n! Warning: \n", color=229)
    for m in msg
        printstyled("    - $m\n", color = 86)
    end
end

"""
Print an error message
"""
function errormsg(msg...)
    !echo_on && return
    printstyled("\n!! ERROR:\n"; color=196)
    for m in msg
        printstyled("    - $m\n", color=86)
    end
end

"""
Print an item.
"""
function item(it)
    !echo_on && return
    printstyled("\n- $it\n"; color=74)
end

"""
Print a `Done` message.
"""
function done(msg... = "Done")
    !echo_on && return
    printstyled(" ... "; color=40)
    printstyled(join(msg, ' '), color=86)
    println()
end

global echo_on = true

set_echo(status::Bool = true) = (global echo_on = status)
