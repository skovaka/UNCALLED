function minimize(name) {
    document.getElementById(name+"-body").style.display = "none";
    document.getElementById(name+"-minimize").style.display = "none";
    document.getElementById(name+"-maximize").style.display = "inline";
}

function maximize(name) {
    document.getElementById(name+"-body").style.display = "block";
    document.getElementById(name+"-minimize").style.display = "inline";
    document.getElementById(name+"-maximize").style.display = "none";
}

function toggle_settings(name) {
    elm = document.getElementById(name+"-settings");
    if (elm.style.display == "none") {
        elm.style.display = "block";
    } else {
        elm.style.display = "none";
    }
}
