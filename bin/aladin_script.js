var cat = A.catalog();
aladin.addCatalog(cat);
aladin.on('objectHovered', function(object) {
    var msg;
    if (object) {
        msg = 'Object ' + object.data.name + ' R.A. = ' + object.ra + ', Dec. = ' + object.dec;
    }
    else {
        msg = ' Move the mouse (or click) on a source for details  ';
    }
    $('#infoDiv').html(msg);
});
// define function triggered when an object is clicked
var objClicked;
aladin.on('objectClicked', function(object) {
    var msg;
    if (object) {
        objClicked = object;
      object.select();
        msg = 'You clicked on object ' + object.data.name + ' located at R.A. = ' + object.ra + ', Dec. = ' + object.dec;
    }
    else {
        objClicked.deselect();
        msg = 'You clicked in void';
    }
    $('#infoDiv').html(msg);
});
