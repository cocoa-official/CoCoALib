-- linkref.lua

function Str(el)
  local s = el.text

  -- [[label|text]]
  local label, text = s:match("%[%[([^|%]]+)|([^%]]+)%]%]")
  if label and text then
    if FORMAT == "html" then
      return pandoc.Link(text, label .. ".html")
    elseif FORMAT == "latex" then
      return pandoc.RawInline("latex", "\\ref{" .. label .. "}")
    else
      return pandoc.Str(text)
    end
  end

  -- [[label]] --> display as code
  label = s:match("%[%[([^%]]+)%]%]")
  if label then
    if FORMAT == "html" then
      return pandoc.Link(pandoc.Code(label), label .. ".html")
    elseif FORMAT == "latex" then
      return pandoc.RawInline("latex", "\\ref{" .. label .. "}")
    else
      return pandoc.Str(label)
    end
  end

  return el
end

