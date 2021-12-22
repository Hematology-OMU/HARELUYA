$(function(){

/**************************************
 * Element 
 *************************************/

/*
 * Disclaimer
 */

var disclaimer = $('<div>');
disclaimer.attr('id', 'disclaimer');
disclaimer.css({
  'position':          'fixed',
  'top':               '50%',
  'left':              '50%',
  'transform':         'translate(-50%,-50%)',
  '-webkit-transform': 'translate(-50%,-50%)',
  'height':            '80%',
  'width':             '80%',
  'z-index':           -1,
  'background-color':  'white',
  'border-color':      'black',
  'border-style':      'solid',
  'border-width':      1,
  'border-radius':     5,
  'text-align':        'center',
  'opacity':           0,
  'padding':           10,
  'overflow-Y': 'scroll'    
});
$('body').append(disclaimer);

var disclaimer_description = $('<img>');
disclaimer_description.attr('id', 'disclaimer_description');
disclaimer_description.attr('src', 'disclaimer.png');
disclaimer_description.css({
  'width':      '100%',
  'max-height': 'auto'
});
disclaimer.append(disclaimer_description);

var accept = $('<button>');
accept.attr('id', 'accept');
accept.text('accept');
accept.css({
  'padding':       10,
  'border-radius': 5
});
$(disclaimer).append(accept);

var cancel = $('<button>');
cancel.attr('id', 'cancel');
cancel.text('cancel');
cancel.css({
  'padding':       10,
  'border-radius': 5
});
$(disclaimer).append(cancel);

/*
 * Curtain
 */

var curtain = $('<div>');
curtain.attr('id', 'curtain');
curtain.css({
  'position':         'absolute',
  'top':              0,
  'left':             0,
  'z-index':          60000,
  'background-color': 'gray',
  'text-align':       'center',
  'opacity':          0.5
});
$('body').append(curtain);

var button = $('<text>');
button.text('Read me');
button.css({
  'position':          'fixed',
  'top':               '50%',
  'left':              '50%',
  'transform':         'translate(-50%,-50%)',
  '-webkit-transform': 'translate(-50%,-50%)',
  'padding':           20,
  'border-radius':     10,
  'font-size':         24,
  'color':             'white',
  'text-shadow':       '0 0 10px white,0 0 15px white',
  'cursor':            'pointer',
  'user-select':       'none'
});
button.hover(function() {
  button.css({
    'color':       'skyblue',
    'text-shadow': '0 0 10px skyblue,0 0 15px skyblue'
  });
}, function() {
  button.css({
    'color':       'white',
    'text-shadow': '0 0 10px white,0 0 15px white'
  });
});
curtain.append(button);

/**************************************
 * Action 
 *************************************/

/*
 * Disclaimer
 */

function showDisclaimer() {
  disclaimer.animate({
    opacity: 1,
    zIndex:  65000
  });
}

function hideDisclaimer() {
  disclaimer.animate({
    opacity: 0,
    zIndex:  -1
  });
}

function toggleDisclaimer() {
  var zIndex = disclaimer.css('z-index');
  if (zIndex === '-1') {
    showDisclaimer();
  } else {
    hideDisclaimer();
  }
}

/*
 * Curtain
 */

function openCurtain() {
  curtain.animate({
    opacity: 0,
    zIndex:  -1
  });
}

/**************************************
 * Bind 
 *************************************/

curtain.on('click', function() {
  toggleDisclaimer();
});

accept.on('click', function() {
  hideDisclaimer();
  openCurtain();
});

cancel.on('click', function() {
  hideDisclaimer();
});

$(window).on('load resize', function() {
  var e = $('body');
  curtain.height(e.height());
  curtain.width(e.width());
});

/**************************************
 * Setup
 *************************************/



});
