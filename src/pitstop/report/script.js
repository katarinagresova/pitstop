document.querySelectorAll("button.fold").forEach(function (b) {
  b.addEventListener("click", function () {
    var body = b.closest("section").querySelector(".sec-b"), open = b.getAttribute("aria-expanded") === "true";
    body.hidden = open; b.setAttribute("aria-expanded", String(!open)); b.textContent = open ? "Show" : "Hide";
  });
});
document.querySelectorAll("button.copy").forEach(function (b) {
  b.addEventListener("click", function () {
    var t = b.parentNode.querySelector("code").textContent;
    if (navigator.clipboard) navigator.clipboard.writeText(t).then(function () {
      b.textContent = "Copied"; setTimeout(function () { b.textContent = "Copy"; }, 1200);
    });
  });
});
