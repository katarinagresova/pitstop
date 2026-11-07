// The site name in the logo's two colours: "pit" in ink, "stop" in the accent (extra.css).
for (const el of document.querySelectorAll(
  ".md-header__topic:first-child .md-ellipsis, .md-nav__title[for='__drawer']")) {
  for (const node of [...el.childNodes]) {
    if (node.nodeType === Node.TEXT_NODE && node.textContent.trim() === "pitstop") {
      const name = document.createElement("span");  // one box, so a flex parent keeps the halves together
      const stop = document.createElement("span");
      stop.className = "rk-stop";
      stop.textContent = "stop";
      name.append("pit", stop);
      node.replaceWith(name);
    }
  }
}
