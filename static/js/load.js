document.getElementById("contact-form").addEventListener("submit", function (event) {
    if (!document.getElementById("input").value.trim()) {
        event.preventDefault();
        var warningMessage = document.getElementById("warning-message");
        warningMessage.style.display = "block";
    } else {
        var warningMessage = document.getElementById("warning-message");
        warningMessage.style.display = "none";
        showLoadingAnimation();
    }
});
