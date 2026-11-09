document.querySelectorAll("button.fold").forEach(function (b) {
  b.addEventListener("click", function () {
    var body = b.closest("section").querySelector(".sec-b"), open = b.getAttribute("aria-expanded") === "true";
    body.hidden = open; b.setAttribute("aria-expanded", String(!open)); b.textContent = open ? "Show" : "Hide";
  });
});
var pick = document.getElementById("pair-pick");
if (pick) {
  var show = function () {
    document.querySelectorAll(".libpanel").forEach(function (p) { p.hidden = p.id !== pick.value; });
  };
  pick.addEventListener("change", show);
  document.querySelectorAll("button[data-step]").forEach(function (b) {
    b.addEventListener("click", function () {
      var i = pick.selectedIndex + Number(b.getAttribute("data-step"));
      if (i >= 0 && i < pick.options.length) { pick.selectedIndex = i; show(); }
    });
  });
  document.querySelectorAll("a[data-pair]").forEach(function (a) {
    a.addEventListener("click", function () { pick.value = a.getAttribute("data-pair"); show(); });
  });
  show();
}
document.querySelectorAll("button.copy").forEach(function (b) {
  b.addEventListener("click", function () {
    var t = b.parentNode.querySelector("code").textContent;
    if (navigator.clipboard) navigator.clipboard.writeText(t).then(function () {
      b.textContent = "Copied"; setTimeout(function () { b.textContent = "Copy"; }, 1200);
    });
  });
});
