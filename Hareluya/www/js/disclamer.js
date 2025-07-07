$(function(){

  // const

  const FRAME_MARGIN = 0.05;

  // member
  
  const screen_ = $('<div>')
    .attr("id", 'screen')
    .css("z-index", "9998")
    .css("position", "absolute")
    .css("top", "0px")
    .css("left", "0px")
    .height(0)
    .width(0)
    .appendTo($('body'));
  
  const disclaimer_ = $('<div>')
    .attr('id', 'disclaimer')
    .css("z-index", "9999")
    .css("position", "absolute")
    .height(0)
    .width(0)
    .appendTo($('body'));

  $('<iframe>')
    .css('width', '100%')
    .css('height', 'calc(100% - 50px)')
    .attr('src', 'disclaimer/disclaimer.html')
    .appendTo(disclaimer_);

  const accept_ = $('<button>')
    .text('accept')
    .appendTo(disclaimer_);

  // action

  function accept() {
    disclaimer_.fadeOut('slow',function () {
      disclaimer_.remove();
    });
    screen_.fadeOut('slow', function () {
      screen_.remove();
    });
    window.onresize = null;
    window.onscroll = null;
  }

  function adjustLayout() {
    const top = window.scrollY;
    const left = window.scrollX;
    const height = window.innerHeight;
    const width = window.innerWidth;
    screen_.offset({ top: top, left: left });
    screen_.height(height);
    screen_.width(width);

    const frame_top = top + (height * FRAME_MARGIN);
    const frame_left = left + (width * FRAME_MARGIN);
    const frame_height = height * (1 - (FRAME_MARGIN * 2));
    const frame_width = width * (1 - (FRAME_MARGIN * 2));
    disclaimer_.offset({ top: frame_top, left: frame_left });
    disclaimer_.height(frame_height);
    disclaimer_.width(frame_width);
  }

  // bind

  accept_.on('click', accept);
  window.onresize = adjustLayout;
  window.onscroll = adjustLayout;

  // init

  adjustLayout();

});
