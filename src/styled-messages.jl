# The color plan below is better in terminal dark mode
"""
Print a title.
"""
function title(msg)
    printstyled('\n', msg, '\n', bold = true, color = :cyan)
    printstyled(repeat('=', length(msg) + 2), '\n', color=31)
end

"""
Print a subtitle.
"""
function subtitle(msg)
    printstyled('\n', msg, '\n', color = :cyan)
    printstyled(repeat('-', length(msg) + 2), '\n', color=31)
end

"""
Print a message.
"""
function message(msg...)
    printstyled("\nMessage: \n", color = :light_magenta)
    for m in msg
        printstyled("    - $m\n", color = 86)
    end
end

"""
Print a warning.
"""
function warning(msg...)
    printstyled("\n! Warning: \n", color=229)
    for m in msg
        printstyled("    - $m\n", color = 86)
    end
end

"""
Print an error message
"""
function errormsg(msg...)
    printstyled("\n!! ERROR:\n"; color=196)
    for m in msg
        printstyled("    - $m\n", color=86)
    end
end

"""
Print an item.
"""
function item(it)
    printstyled("\n- $it\n"; color=74)
end

"""
Print a `Done` message.
"""
function done(msg... = "Done")
    printstyled(" ... "; color=40)
    printstyled(join(msg, ' '), color=86)
    println()
end
