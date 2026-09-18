"use strict";
const select = document.querySelector("#side-vertices");
const cards = [...document.querySelectorAll(".side-card")];
select.addEventListener("change", () => {
  for (const card of cards) card.hidden = select.value !== "all" && card.dataset.sideVertices !== select.value;
  for (const group of document.querySelectorAll(".side-group")) {
    group.hidden = [...group.querySelectorAll(".side-card")].every(card => card.hidden);
    document.querySelector(`.side-toolbar a[href="#${group.id}"]`).hidden = group.hidden;
  }
  const shown = cards.filter(card => !card.hidden).length;
  document.querySelector("#side-status").textContent = `Showing ${shown} of ${cards.length} patterns, ordered by quad count, then best quality first.`;
});
